###############################################################################
## cengiz_application.R
## Apply the l0 (BIC) and DP-Bayes specification estimators to the event-level
## treatment effects from the stacked minimum-wage design of Cengiz et al. (2019).
##
## INPUT (summary-statistic form -- preferred, minimal, and EXACT here):
##   a CSV `cengiz_events.csv` with one row per event and columns
##     event : event/stack identifier (138 rows)
##     tau   : flexible (per-event) CATT estimate  tau_hat_g
##     se    : its standard error (clustered as in the original study)
##   Optionally `weight` (e.g. affected employment) for the overall ATT.
##
## Why summary statistics suffice: the stacked design is block-diagonal across
## events (each event uses its own treated state and clean controls over its own
## window), so the 138 CATT estimates are mutually (near-)independent and
## tau_hat_g ~ N(tau_g, se_g^2).  The l0/DP clustering therefore operates
## exactly on (tau_hat_g, se_g) -- no cross-event covariance to model, in
## contrast to the fixed-effects-correlated simulation.
##
## If you instead have the stacked micro panel, produce cengiz_events.csv first
## by estimating the fully flexible (event-interacted) TWFE -- see the commented
## fixest template at the bottom -- then run this script.
###############################################################################

infile  <- "cengiz_events.csv"     # <-- point this at your event-level file
outstub <- "cengiz"

if (!file.exists(infile)) {
  stop(sprintf(paste0(
    "Input '%s' not found.\n",
    "Provide a CSV with columns event, tau, se (one row per event; ~138 rows).\n",
    "If you only have the stacked micro panel, first estimate the flexible\n",
    "event-interacted TWFE to obtain per-event CATTs + SEs (see template below)."), infile))
}

ev <- read.csv(infile, stringsAsFactors = FALSE)
stopifnot(all(c("tau","se") %in% names(ev)))
if (is.null(ev$event)) ev$event <- seq_len(nrow(ev))
if (is.null(ev$weight)) ev$weight <- 1

## --- data hygiene: flag events with degenerate normalization -------------
## Cengiz et al. normalize by the pre-period share of affected employment;
## very thin cells give a near-zero denominator and explosive estimates
## (the paper flags three FIPS-1 events). Winsorize |tau| at a robust cutoff.
mad_tau <- mad(ev$tau); med_tau <- median(ev$tau)
cutoff  <- med_tau + 8 * mad_tau          # generous; adjust as needed
flagged <- which(abs(ev$tau - med_tau) > 8 * mad_tau | ev$se <= 0 | !is.finite(ev$tau))
if (length(flagged)) {
  cat(sprintf("Flagged %d event(s) with anomalous magnitude/SE (dropped): %s\n",
              length(flagged), paste(ev$event[flagged], collapse = ", ")))
  ev <- ev[-flagged, , drop = FALSE]
}
E    <- nrow(ev)
tauh <- ev$tau
prec <- 1 / ev$se^2                       # precision weights n_g = 1/se_g^2
wt   <- ev$weight / sum(ev$weight)
cat(sprintf("Using %d events.\n", E))

## precision-weighted ("pooled") and simple flexible ATT summaries
att_pool <- sum(prec * tauh) / sum(prec)  # single common effect (pooled TWFE)
att_flex <- sum(wt * tauh)                # weighted average of flexible CATTs

## ------------------------------- helpers ----------------------------------
memb <- function(labels) {
  u <- sort(unique(labels)); A <- matrix(0, length(labels), length(u))
  for (j in seq_along(u)) A[labels == u[j], j] <- 1; A
}
## grouped weighted RSS and effects given a label vector
grp_fit <- function(labels) {
  phi <- tapply(prec * tauh, labels, sum) / tapply(prec, labels, sum)
  tau_g <- phi[as.character(labels)]
  rss <- sum(prec * (tauh - tau_g)^2)     # weighted (chi-square-scaled) deviance
  list(tau = as.numeric(tau_g), phi = phi, rss = rss, m = length(phi))
}

## ---- l0: greedy agglomeration path + BIC (reduced-model, n = E) ----------
l0_events <- function() {
  groups <- as.list(seq_len(E)); g_tau <- tauh; g_n <- prec
  labs_at <- vector("list", E); labs_at[[E]] <- seq_len(E)
  for (step in seq_len(E - 1)) {
    ng <- length(groups); best <- Inf; bi <- bj <- NA
    for (a in 1:(ng - 1)) for (b in (a + 1):ng) {
      cost <- (g_n[a]*g_n[b]/(g_n[a]+g_n[b]))*(g_tau[a]-g_tau[b])^2
      if (cost < best) { best <- cost; bi <- a; bj <- b }
    }
    nn <- g_n[bi]+g_n[bj]; tt <- (g_n[bi]*g_tau[bi]+g_n[bj]*g_tau[bj])/nn
    groups[[bi]] <- c(groups[[bi]], groups[[bj]]); g_tau[bi] <- tt; g_n[bi] <- nn
    groups[[bj]] <- NULL; g_tau <- g_tau[-bj]; g_n <- g_n[-bj]
    labs <- integer(E); for (gi in seq_along(groups)) labs[groups[[gi]]] <- gi
    labs_at[[E - step]] <- labs
  }
  bic <- sapply(1:E, function(m) { f <- grp_fit(labs_at[[m]]); f$rss + m*log(E) })
  ms  <- which.min(bic)
  list(labels = labs_at[[ms]], m = ms, fit = grp_fit(labs_at[[ms]]))
}

## ---- DP normal mixture with known variances (collapsed Gibbs) ------------
dp_events <- function(alpha = max(1, E/8), mu0 = median(tauh),
                      s0sq = (4*mad(tauh))^2, iters = 2000, burn = 500) {
  z <- rep(1L, E); b0 <- 1/s0sq
  post_tau <- numeric(E); coclust <- matrix(0, E, E); kept <- 0
  att_dr <- numeric(iters - burn); m_dr <- numeric(iters - burn)
  for (it in seq_len(iters)) {
    for (k in seq_len(E)) {
      z[k] <- NA; cl <- sort(unique(z[!is.na(z)])); logp <- numeric(length(cl)+1)
      for (j in seq_along(cl)) { mem <- which(z==cl[j]); Pc <- b0+sum(prec[mem])
        mcj <- (mu0*b0+sum(prec[mem]*tauh[mem]))/Pc
        logp[j] <- log(length(mem)) + dnorm(tauh[k], mcj, sqrt(1/prec[k]+1/Pc), log=TRUE) }
      logp[length(cl)+1] <- log(alpha) + dnorm(tauh[k], mu0, sqrt(1/prec[k]+s0sq), log=TRUE)
      p <- exp(logp-max(logp)); p <- p/sum(p); pk <- sample.int(length(p),1,prob=p)
      z[k] <- if (pk<=length(cl)) cl[pk] else max(cl,0)+1
    }
    z <- match(z, sort(unique(z)))
    if (it > burn) { kept <- kept+1; phi <- numeric(max(z))
      for (j in seq_len(max(z))) { mem <- which(z==j); Pc <- b0+sum(prec[mem])
        mcj <- (mu0*b0+sum(prec[mem]*tauh[mem]))/Pc; phi[j] <- rnorm(1, mcj, sqrt(1/Pc)) }
      post_tau <- post_tau + phi[z]; coclust <- coclust + outer(z, z, `==`)
      att_dr[kept] <- sum(wt * phi[z]); m_dr[kept] <- max(z) }
  }
  list(tau = post_tau/kept, coclust = coclust/kept, labels = z,
       m_mean = mean(m_dr), att = mean(att_dr),
       att_lo = quantile(att_dr, .025), att_hi = quantile(att_dr, .975))
}

## --------------------------------- run ------------------------------------
l0 <- l0_events()
dp <- dp_events()

cat("\n================= Cengiz application: specification =================\n")
cat(sprintf("Events (after hygiene): %d\n", E))
cat(sprintf("Pooled TWFE effect (precision-weighted common effect): %.4f\n", att_pool))
cat(sprintf("Flexible ATT (weighted avg of %d CATTs):               %.4f\n", E, att_flex))
cat(sprintf("l0 (BIC):  %d groups; overall ATT = %.4f\n", l0$m, sum(wt*l0$fit$tau)))
cat(sprintf("DP Bayes:  %.1f groups (posterior mean); ATT = %.4f  [%.4f, %.4f]\n",
            dp$m_mean, dp$att, dp$att_lo, dp$att_hi))

## save per-event assignments and grouped effects
res <- data.frame(event = ev$event, tau = tauh, se = ev$se,
                  l0_group = l0$labels, l0_effect = l0$fit$tau,
                  dp_group = dp$labels, dp_effect = dp$tau)
write.csv(res, paste0(outstub, "_events_grouped.csv"), row.names = FALSE)

## co-clustering (posterior similarity) heatmap
png(paste0(outstub, "_coclustering.png"), 1400, 1300, res = 200)
ord <- order(dp$tau)
image(dp$coclust[ord, ord], axes = FALSE,
      main = "DP posterior co-clustering probability (events ordered by effect)",
      col = grey.colors(20, start = 1, end = 0))
dev.off()

## effect-by-event with l0 groups
png(paste0(outstub, "_effects.png"), 1700, 900, res = 200)
o <- order(tauh)
plot(seq_len(E), tauh[o], pch = 19, col = l0$labels[o] + 1,
     xlab = "event (sorted by flexible CATT)", ylab = "treatment effect",
     main = sprintf("Flexible CATTs and l0 grouping (%d groups)", l0$m))
arrows(seq_len(E), (tauh - 1.96*ev$se)[o], seq_len(E), (tauh + 1.96*ev$se)[o],
       length = 0, col = "grey70")
points(seq_len(E), l0$fit$tau[o], pch = 4, col = "blue")
abline(h = att_pool, lty = 2)
legend("topleft", c("flexible CATT", "l0 group effect", "pooled"),
       pch = c(19, 4, NA), lty = c(NA, NA, 2), col = c("black","blue","black"), bty = "n")
dev.off()

cat(sprintf("\nWrote %s_events_grouped.csv, %s_coclustering.png, %s_effects.png\n",
            outstub, outstub, outstub))

###############################################################################
## TEMPLATE: building cengiz_events.csv from the stacked micro panel (fixest)
## -- uncomment and adapt column names to your stacked data `dat`:
##
## library(fixest); library(broom)
## # dat columns: y (outcome), treatpost (D_itg), state, period, stack (event g)
## m_flex <- feols(y ~ i(stack, treatpost, ref = 0) | state^stack + period^stack,
##                 data = dat, cluster = ~ state)
## ct <- broom::tidy(m_flex)
## ct <- ct[grepl("treatpost", ct$term), ]
## ev <- data.frame(event = seq_len(nrow(ct)), tau = ct$estimate, se = ct$std.error)
## write.csv(ev, "cengiz_events.csv", row.names = FALSE)
###############################################################################
