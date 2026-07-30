###############################################################################
## staggered_sim.R
## Monte Carlo for "Specification in TWFE in Staggered Timing and Partial
## Heterogeneous Effect".  Variance-first evaluation.
##
## PART 1  delta-sweep x m-grid : sampling variance / bias / recovery
## PART 2  coverage (E6)        : 95% CI/credible-interval coverage & length,
##                               incl. naive vs sample-split l0 and DP intervals
## PART 3  solution paths        : ATT bias/variance and CATT variance as a
##                               function of the l0 penalty (via m on the greedy
##                               path) and of the DP concentration alpha
##
## Design: balanced, orthogonal (one cell per cohort-time). Linear estimators
## reduce to operations on the K x K within-Gram matrix G.
###############################################################################

dir.create("figures", showWarnings = FALSE)
set.seed(1234)

## ----------------------------- configuration ------------------------------
N        <- 2000
T        <- 10
cohorts  <- c(0, 3, 5, 7)
m_grid   <- c(1, 3, 6, 9, 18)
delta_grid <- c(3, 6, 12)        # NEW: separation sweep
sigma    <- 1

R_reps_main <- 500               # PART 1
R_reps_cov  <- 500               # PART 2
R_reps_path <- 500               # PART 3

## DP-Bayes hyperparameters
dp_mu0   <- 0
dp_s0sq  <- 100
dp_iters <- 150
dp_burn  <- 50
a0_sig   <- 2                    # Inv-Gamma(a0,b0) prior on sigma^2 (full Bayesian)
b0_sig   <- 1
dp_alpha_default <- 7            # E[m] ~ K/2
alpha_grid <- c(1, 3, 7, 15, 30) # PART 3 DP path

out_dir  <- "."

## --------------------------- shared structure -----------------------------
stopifnot(N %% length(cohorts) == 0)
cohort_full <- rep(cohorts, each = N / length(cohorts))
cells <- do.call(rbind, lapply(cohorts[cohorts > 0], function(g)
  data.frame(g = g, t = g:T)))
K <- nrow(cells)
w_att <- rep(1 / K, K)           # equal-weight overall ATT aggregator
cat(sprintf("Design: N=%d, T=%d, K=%d cells\n", N, T, K))

demean2 <- function(M) M - rowMeans(M) - rep(colMeans(M), each = nrow(M)) + mean(M)

## build within-transformed design for a given per-unit cohort vector
build_design <- function(coh) {
  n <- length(coh); NT <- n * T
  Dt <- matrix(0, NT, K)
  for (k in seq_len(K)) {
    Dk <- matrix(0, n, T); Dk[coh == cells$g[k], cells$t[k]] <- 1
    Dt[, k] <- as.vector(demean2(Dk))
  }
  G <- crossprod(Dt)
  list(Dt = Dt, G = G, Ginv = solve(G), nk = diag(G), n = n, NT = NT)
}

des_full <- build_design(cohort_full)
sd_flex_bar <- mean(sigma / sqrt(des_full$nk))
cat(sprintf("Mean flexible sd = %.4f\n", sd_flex_bar))

## fixed 50/50 unit split (stratified by cohort) for sample-split inference
idxA <- unlist(lapply(cohorts, function(g) {
  ii <- which(cohort_full == g); ii[seq_len(length(ii) / 2)]
}))
splitA <- rep(FALSE, N); splitA[idxA] <- TRUE
des_A <- build_design(cohort_full[splitA])
des_B <- build_design(cohort_full[!splitA])

## ------------------------------- helpers ----------------------------------
make_truth <- function(m, Delta) {
  sizes <- rep(K %/% m, m); rem <- K %% m
  if (rem > 0) sizes[seq_len(rem)] <- sizes[seq_len(rem)] + 1
  labels <- rep(seq_len(m), times = sizes)
  mus <- (seq_len(m) - (m + 1) / 2) * Delta
  list(labels = labels, tau = mus[labels])
}
memb <- function(labels) {
  u <- sort(unique(labels)); A <- matrix(0, length(labels), length(u))
  for (j in seq_along(u)) A[labels == u[j], j] <- 1
  A
}
fit_group <- function(labels, c_vec, yty, Gmat) {
  A <- memb(labels); GA <- crossprod(A, Gmat %*% A); cA <- crossprod(A, c_vec)
  GAi <- solve(GA); phi <- GAi %*% cA
  list(tau = as.vector(A %*% phi), rss = yty - sum(cA * phi),
       m = ncol(A), A = A, GAi = GAi)
}
## per-cell standard errors for a grouping, given design and variance
cell_se <- function(fit, Gmat, sigma2) {
  V <- sigma2 * (fit$A %*% fit$GAi %*% t(fit$A))
  sqrt(pmax(diag(V), 0))
}
att_se <- function(fit, sigma2) {
  V <- sigma2 * (fit$A %*% fit$GAi %*% t(fit$A))
  sqrt(as.numeric(t(w_att) %*% V %*% w_att))
}
sig2_hat <- function(rss, des, m) rss / (des$NT - (des$n - 1) - (T - 1) - m)

## full l0 greedy path + BIC
l0_run <- function(des, c_vec, yty, tau_flex) {
  groups <- as.list(seq_len(K)); g_tau <- tau_flex; g_n <- des$nk
  labels_at <- vector("list", K); labels_at[[K]] <- seq_len(K)
  for (step in seq_len(K - 1)) {
    ng <- length(groups); best <- Inf; bi <- bj <- NA
    for (a in 1:(ng - 1)) for (b in (a + 1):ng) {
      cost <- (g_n[a] * g_n[b] / (g_n[a] + g_n[b])) * (g_tau[a] - g_tau[b])^2
      if (cost < best) { best <- cost; bi <- a; bj <- b }
    }
    merged <- c(groups[[bi]], groups[[bj]]); new_n <- g_n[bi] + g_n[bj]
    new_t <- (g_n[bi] * g_tau[bi] + g_n[bj] * g_tau[bj]) / new_n
    groups[[bi]] <- merged; g_tau[bi] <- new_t; g_n[bi] <- new_n
    groups[[bj]] <- NULL; g_tau <- g_tau[-bj]; g_n <- g_n[-bj]
    labs <- integer(K); for (gi in seq_along(groups)) labs[groups[[gi]]] <- gi
    labels_at[[K - step]] <- labs
  }
  fits <- vector("list", K); bic <- numeric(K)
  for (m in 1:K) {
    f <- fit_group(labels_at[[m]], c_vec, yty, des$G); fits[[m]] <- f
    bic[m] <- des$NT * log(f$rss / des$NT) + m * log(des$NT)
  }
  list(labels_at = labels_at, fits = fits, bic = bic, m_sel = which.min(bic))
}

## DP collapsed Gibbs on cell sufficient stats; returns mean, labels, per-cell
## and ATT credible-interval quantiles
## Label moves use the fast diagonal (independence) approximation for the
## marginal likelihood; the RECORDED draws use the EXACT GLS posterior with the
## full within covariance A'GA (paper eqs 40-41), so credible intervals reflect
## the true cross-cell correlation induced by the fixed effects.
dp_bayes <- function(tau_flex, c_vec, des, alpha = dp_alpha_default) {
  prec <- des$nk / sigma^2
  z <- rep(1L, K); b0 <- 1 / dp_s0sq
  draws <- matrix(NA, dp_iters - dp_burn, K); att_dr <- numeric(dp_iters - dp_burn)
  kept <- 0
  for (it in seq_len(dp_iters)) {
    for (k in seq_len(K)) {
      z[k] <- NA; cl <- sort(unique(z[!is.na(z)]))
      logp <- numeric(length(cl) + 1)
      for (j in seq_along(cl)) {
        mem <- which(z == cl[j]); Pc <- b0 + sum(prec[mem])
        mcj <- (dp_mu0 * b0 + sum(prec[mem] * tau_flex[mem])) / Pc
        logp[j] <- log(length(mem)) +
          dnorm(tau_flex[k], mcj, sqrt(1 / prec[k] + 1 / Pc), log = TRUE)
      }
      logp[length(cl) + 1] <- log(alpha) +
        dnorm(tau_flex[k], dp_mu0, sqrt(1 / prec[k] + dp_s0sq), log = TRUE)
      p <- exp(logp - max(logp)); p <- p / sum(p)
      pick <- sample.int(length(p), 1, prob = p)
      z[k] <- if (pick <= length(cl)) cl[pick] else (max(cl, 0) + 1)
    }
    z <- match(z, sort(unique(z)))
    if (it > dp_burn) {
      kept <- kept + 1
      A <- memb(z); mm <- ncol(A)
      GA <- crossprod(A, des$G %*% A); cA <- crossprod(A, c_vec)
      Vs <- solve(GA / sigma^2 + diag(mm) / dp_s0sq)
      mu <- Vs %*% (cA / sigma^2 + dp_mu0 / dp_s0sq)
      phi <- as.vector(mu + t(chol(Vs)) %*% rnorm(mm))
      tau <- as.vector(A %*% phi)
      draws[kept, ] <- tau; att_dr[kept] <- sum(w_att * tau)
    }
  }
  qc <- apply(draws, 2, quantile, c(0.025, 0.975))
  list(tau = colMeans(draws), labels = z, m = length(unique(z)),
       lo = qc[1, ], hi = qc[2, ],
       att = mean(att_dr), att_lo = quantile(att_dr, 0.025),
       att_hi = quantile(att_dr, 0.975))
}

## FULL Bayesian DP: sigma^2 given an inverse-gamma prior and drawn each sweep
## from its conjugate full conditional (residual df), then cell precisions are
## updated.  Recorded phi uses the EXACT GLS posterior with A'GA, as above.  This
## is the PRIMARY estimator; the fixed-sigma^2 dp_bayes() above is the special
## case used only for the appendix comparison and the exact l0 correspondence.
dp_bayes_rand <- function(tau_flex, c_vec, yty, des, alpha = dp_alpha_default,
                          a0 = a0_sig, b0 = b0_sig) {
  z <- rep(1L, K); b0m <- 1 / dp_s0sq; sig2 <- 1
  nkeep <- dp_iters - dp_burn
  draws <- matrix(NA, nkeep, K); att_dr <- numeric(nkeep); sig_dr <- numeric(nkeep)
  kept <- 0
  for (it in seq_len(dp_iters)) {
    prec <- des$nk / sig2
    for (k in seq_len(K)) {
      z[k] <- NA; cl <- sort(unique(z[!is.na(z)]))
      logp <- numeric(length(cl) + 1)
      for (j in seq_along(cl)) {
        mem <- which(z == cl[j]); Pc <- b0m + sum(prec[mem])
        mcj <- (dp_mu0 * b0m + sum(prec[mem] * tau_flex[mem])) / Pc
        logp[j] <- log(length(mem)) +
          dnorm(tau_flex[k], mcj, sqrt(1 / prec[k] + 1 / Pc), log = TRUE)
      }
      logp[length(cl) + 1] <- log(alpha) +
        dnorm(tau_flex[k], dp_mu0, sqrt(1 / prec[k] + dp_s0sq), log = TRUE)
      p <- exp(logp - max(logp)); p <- p / sum(p)
      pick <- sample.int(length(p), 1, prob = p)
      z[k] <- if (pick <= length(cl)) cl[pick] else (max(cl, 0) + 1)
    }
    z <- match(z, sort(unique(z)))
    ## exact GLS posterior for phi given current sig2
    A  <- memb(z); mm <- ncol(A)
    GA <- crossprod(A, des$G %*% A); cA <- crossprod(A, c_vec)
    Vs <- solve(GA / sig2 + diag(mm) / dp_s0sq)
    mu <- Vs %*% (cA / sig2 + dp_mu0 / dp_s0sq)
    phi <- as.vector(mu + t(chol(Vs)) %*% rnorm(mm))
    tau <- as.vector(A %*% phi)
    ## draw sigma^2 from its inverse-gamma full conditional (residual df)
    rss_full <- yty - 2 * sum(cA * phi) + as.numeric(t(phi) %*% GA %*% phi)
    rss_full <- max(rss_full, 1e-8)
    df_res <- des$NT - (des$n - 1) - (T - 1) - mm
    sig2 <- 1 / rgamma(1, a0 + df_res / 2, b0 + rss_full / 2)
    if (it > dp_burn) {
      kept <- kept + 1
      draws[kept, ] <- tau; att_dr[kept] <- sum(w_att * tau); sig_dr[kept] <- sig2
    }
  }
  qc <- apply(draws, 2, quantile, c(0.025, 0.975))
  list(tau = colMeans(draws), labels = z, m = length(unique(z)),
       lo = qc[1, ], hi = qc[2, ], sig2 = mean(sig_dr),
       att = mean(att_dr), att_lo = quantile(att_dr, 0.025),
       att_hi = quantile(att_dr, 0.975))
}

adj_rand <- function(a, b) {
  tab <- table(a, b); n <- sum(tab); comb2 <- function(x) x * (x - 1) / 2
  si <- sum(comb2(tab)); sa <- sum(comb2(rowSums(tab))); sb <- sum(comb2(colSums(tab)))
  ex <- sa * sb / comb2(n); mx <- (sa + sb) / 2
  if (mx == ex) return(1); (si - ex) / (mx - ex)
}

## draw one within-transformed replication for a design and truth
draw_rep <- function(des, coh, tau_true) {
  TAU <- matrix(0, des$n, T)
  for (k in seq_len(K)) TAU[coh == cells$g[k], cells$t[k]] <- tau_true[k]
  ytil <- as.vector(demean2(TAU + matrix(rnorm(des$n * T, 0, sigma), des$n, T)))
  list(c = crossprod(des$Dt, ytil), yty = sum(ytil^2))
}

has_genlasso <- requireNamespace("genlasso", quietly = TRUE)
cat(sprintf("genlasso available: %s\n", has_genlasso))

###############################################################################
## PART 1 : delta-sweep x m-grid  (variance-first)
###############################################################################
cat("\n===== PART 1: delta-sweep x m-grid =====\n")
methods1 <- c("Pooled", "Flexible", "Oracle_PH", "l0_BIC", "DP_Bayes")
res1 <- list()
for (dt in delta_grid) {
  Delta <- dt * sd_flex_bar
  for (m_true in m_grid) {
    truth <- make_truth(m_true, Delta); tau_true <- truth$tau; lab_true <- truth$labels
    est <- lapply(methods1, function(x) matrix(NA, R_reps_main, K)); names(est) <- methods1
    mh <- ar <- matrix(NA, R_reps_main, length(methods1)); colnames(mh) <- colnames(ar) <- methods1
    for (r in seq_len(R_reps_main)) {
      d <- draw_rep(des_full, cohort_full, tau_true); c_vec <- d$c; yty <- d$yty
      tf <- as.vector(des_full$Ginv %*% c_vec)
      est$Flexible[r, ]  <- tf
      est$Pooled[r, ]    <- fit_group(rep(1L, K), c_vec, yty, des_full$G)$tau
      est$Oracle_PH[r, ] <- fit_group(lab_true, c_vec, yty, des_full$G)$tau
      mh[r, "Flexible"] <- K; mh[r, "Pooled"] <- 1; mh[r, "Oracle_PH"] <- m_true
      ar[r, "Flexible"] <- adj_rand(seq_len(K), lab_true)
      ar[r, "Pooled"]   <- adj_rand(rep(1, K), lab_true); ar[r, "Oracle_PH"] <- 1
      l0 <- l0_run(des_full, c_vec, yty, tf); ms <- l0$m_sel
      est$l0_BIC[r, ] <- l0$fits[[ms]]$tau
      mh[r, "l0_BIC"] <- ms; ar[r, "l0_BIC"] <- adj_rand(l0$labels_at[[ms]], lab_true)
      dp <- dp_bayes_rand(tf, c_vec, yty, des_full)
      est$DP_Bayes[r, ] <- dp$tau; mh[r, "DP_Bayes"] <- dp$m
      ar[r, "DP_Bayes"] <- adj_rand(dp$labels, lab_true)
    }
    for (mth in methods1) {
      E <- est[[mth]]; cv <- apply(E, 2, var); cb <- colMeans(E) - tau_true
      res1[[length(res1) + 1]] <- data.frame(
        delta = dt, m_true = m_true, method = mth,
        variance = mean(cv), bias_abs = mean(abs(cb)),
        rmse = mean(sqrt(cb^2 + cv)), mhat = mean(mh[, mth]), ari = mean(ar[, mth]))
    }
    cat(sprintf("  delta=%2d  m*=%2d done\n", dt, m_true))
  }
}
res1 <- do.call(rbind, res1)
for (dt in delta_grid) for (m_true in m_grid) {
  vf <- res1$variance[res1$delta == dt & res1$m_true == m_true & res1$method == "Flexible"]
  sel <- res1$delta == dt & res1$m_true == m_true
  res1$var_ratio[sel] <- res1$variance[sel] / vf
}
res1[, c("variance","bias_abs","rmse")] <- round(res1[, c("variance","bias_abs","rmse")], 5)
res1[, c("mhat","ari","var_ratio")] <- round(res1[, c("mhat","ari","var_ratio")], 3)
write.csv(res1, file.path(out_dir, "sim_results.csv"), row.names = FALSE)

## Figure: variance vs m* at headline delta=6
png(file.path("figures", "variance_vs_mstar.png"), 1600, 1100, res = 200)
sub6 <- res1[res1$delta == 6, ]; cols <- setNames(seq_along(methods1) + 1, methods1)
plot(NA, xlim = range(m_grid), ylim = range(sub6$variance),
     xlab = expression(m^"*"), ylab = "Avg sampling variance of CATT",
     main = expression("Variance under partial homogeneity ("*delta*"=6)"))
for (mth in methods1) { s <- sub6[sub6$method == mth, ]
  lines(s$m_true, s$variance, col = cols[mth], lwd = 2, type = "b", pch = 19) }
legend("topleft", methods1, col = cols, lwd = 2, pch = 19, bty = "n")
dev.off()

## Figure: variance ratio & ARI vs delta at m*=6 (separation dependence)
png(file.path("figures", "separation_vs_delta.png"), 1800, 900, res = 200)
par(mfrow = c(1, 2)); fm <- c("Oracle_PH","l0_BIC","DP_Bayes")
s6 <- res1[res1$m_true == 6 & res1$method %in% fm, ]
plot(NA, xlim = range(delta_grid), ylim = c(0, max(s6$var_ratio)),
     xlab = expression(delta), ylab = "Variance ratio to Flexible", main = "m*=6: variance")
abline(h = 1, lty = 3)
for (mth in fm) { s <- s6[s6$method == mth, ]; lines(s$delta, s$var_ratio, col = cols[mth], lwd = 2, type = "b", pch = 19) }
legend("topright", fm, col = cols[fm], lwd = 2, pch = 19, bty = "n")
plot(NA, xlim = range(delta_grid), ylim = c(0, 1),
     xlab = expression(delta), ylab = "Adjusted Rand Index", main = "m*=6: recovery")
for (mth in c("l0_BIC","DP_Bayes")) { s <- s6[s6$method == mth, ]; lines(s$delta, s$ari, col = cols[mth], lwd = 2, type = "b", pch = 19) }
legend("bottomright", c("l0_BIC","DP_Bayes"), col = cols[c("l0_BIC","DP_Bayes")], lwd = 2, pch = 19, bty = "n")
dev.off()

cat("PART 1 done -> sim_results.csv, variance_vs_mstar.png, separation_vs_delta.png\n")

###############################################################################
## PART 2 : coverage experiment (E6) at m*=6, delta=6
###############################################################################
cat("\n===== PART 2: coverage (m*=6, delta=6) =====\n")
Delta <- 6 * sd_flex_bar
truth <- make_truth(6, Delta); tau_true <- truth$tau; att_true <- sum(w_att * tau_true)
methods2 <- c("Pooled","Flexible","Oracle_PH","l0_naive","l0_split","DP_Bayes","DP_fixed")
cov_cat <- cov_att <- len_cat <- len_att <- matrix(0, R_reps_cov, length(methods2))
colnames(cov_cat) <- colnames(cov_att) <- colnames(len_cat) <- colnames(len_att) <- methods2
z95 <- qnorm(0.975)
ci_stats <- function(tau_hat, se, att_hat, att_se_v) {
  list(cc = mean(tau_true >= tau_hat - z95*se & tau_true <= tau_hat + z95*se),
       lc = mean(2*z95*se),
       ca = as.numeric(att_true >= att_hat - z95*att_se_v & att_true <= att_hat + z95*att_se_v),
       la = 2*z95*att_se_v)
}
for (r in seq_len(R_reps_cov)) {
  d <- draw_rep(des_full, cohort_full, tau_true); c_vec <- d$c; yty <- d$yty
  tf <- as.vector(des_full$Ginv %*% c_vec)
  ## Flexible
  ff <- fit_group(seq_len(K), c_vec, yty, des_full$G); s2 <- sig2_hat(ff$rss, des_full, K)
  se <- cell_se(ff, des_full$G, s2); st <- ci_stats(tf, se, sum(w_att*tf), att_se(ff, s2))
  cov_cat[r,"Flexible"]<-st$cc; len_cat[r,"Flexible"]<-st$lc; cov_att[r,"Flexible"]<-st$ca; len_att[r,"Flexible"]<-st$la
  ## Pooled
  fp <- fit_group(rep(1L,K), c_vec, yty, des_full$G); s2p <- sig2_hat(fp$rss, des_full, 1)
  se <- cell_se(fp, des_full$G, s2p); st <- ci_stats(fp$tau, se, sum(w_att*fp$tau), att_se(fp, s2p))
  cov_cat[r,"Pooled"]<-st$cc; len_cat[r,"Pooled"]<-st$lc; cov_att[r,"Pooled"]<-st$ca; len_att[r,"Pooled"]<-st$la
  ## Oracle
  fo <- fit_group(truth$labels, c_vec, yty, des_full$G); s2o <- sig2_hat(fo$rss, des_full, 6)
  se <- cell_se(fo, des_full$G, s2o); st <- ci_stats(fo$tau, se, sum(w_att*fo$tau), att_se(fo, s2o))
  cov_cat[r,"Oracle_PH"]<-st$cc; len_cat[r,"Oracle_PH"]<-st$lc; cov_att[r,"Oracle_PH"]<-st$ca; len_att[r,"Oracle_PH"]<-st$la
  ## l0 naive: select + estimate + CI all on full data (partition treated fixed)
  l0 <- l0_run(des_full, c_vec, yty, tf); ms <- l0$m_sel; fl <- l0$fits[[ms]]
  s2l <- sig2_hat(fl$rss, des_full, ms); se <- cell_se(fl, des_full$G, s2l)
  st <- ci_stats(fl$tau, se, sum(w_att*fl$tau), att_se(fl, s2l))
  cov_cat[r,"l0_naive"]<-st$cc; len_cat[r,"l0_naive"]<-st$lc; cov_att[r,"l0_naive"]<-st$ca; len_att[r,"l0_naive"]<-st$la
  ## l0 sample-split: select partition on A, estimate + CI on B
  dA <- draw_rep(des_A, cohort_full[splitA], tau_true)
  tfA <- as.vector(des_A$Ginv %*% dA$c)
  l0A <- l0_run(des_A, dA$c, dA$yty, tfA); labA <- l0A$labels_at[[l0A$m_sel]]
  dB <- draw_rep(des_B, cohort_full[!splitA], tau_true)
  fB <- fit_group(labA, dB$c, dB$yty, des_B$G); s2B <- sig2_hat(fB$rss, des_B, length(unique(labA)))
  se <- cell_se(fB, des_B$G, s2B); st <- ci_stats(fB$tau, se, sum(w_att*fB$tau), att_se(fB, s2B))
  cov_cat[r,"l0_split"]<-st$cc; len_cat[r,"l0_split"]<-st$lc; cov_att[r,"l0_split"]<-st$ca; len_att[r,"l0_split"]<-st$la
  ## DP Bayes credible intervals: full Bayesian (random sigma^2) = PRIMARY
  dp <- dp_bayes_rand(tf, c_vec, yty, des_full)
  cov_cat[r,"DP_Bayes"]<-mean(tau_true>=dp$lo & tau_true<=dp$hi); len_cat[r,"DP_Bayes"]<-mean(dp$hi-dp$lo)
  cov_att[r,"DP_Bayes"]<-as.numeric(att_true>=dp$att_lo & att_true<=dp$att_hi); len_att[r,"DP_Bayes"]<-dp$att_hi-dp$att_lo
  ## DP with fixed (plug-in) sigma^2 = APPENDIX comparison
  dpf <- dp_bayes(tf, c_vec, des_full)
  cov_cat[r,"DP_fixed"]<-mean(tau_true>=dpf$lo & tau_true<=dpf$hi); len_cat[r,"DP_fixed"]<-mean(dpf$hi-dpf$lo)
  cov_att[r,"DP_fixed"]<-as.numeric(att_true>=dpf$att_lo & att_true<=dpf$att_hi); len_att[r,"DP_fixed"]<-dpf$att_hi-dpf$att_lo
}
cov_tab <- data.frame(method = methods2,
  cover_CATT = round(colMeans(cov_cat), 3), len_CATT = round(colMeans(len_cat), 4),
  cover_ATT = round(colMeans(cov_att), 3), len_ATT = round(colMeans(len_att), 4))
write.csv(cov_tab, file.path(out_dir, "coverage_results.csv"), row.names = FALSE)
cat("Nominal 95% coverage (target = true CATTs / true ATT):\n"); print(cov_tab, row.names = FALSE)

###############################################################################
## PART 3 : solution paths at m*=6, delta=6
###############################################################################
cat("\n===== PART 3: solution paths (m*=6, delta=6) =====\n")
## l0 path: index by m on the greedy path (m large = small penalty)
l0_att_b <- l0_att_v <- l0_catt_v <- matrix(NA, R_reps_path, K)
msel_l0 <- numeric(R_reps_path)
## DP path: index by alpha
dp_att_b <- dp_att_v <- dp_catt_v <- matrix(NA, R_reps_path, length(alpha_grid))
dp_att_store <- array(NA, c(R_reps_path, length(alpha_grid)))
dp_att_mat <- array(NA, c(R_reps_path, length(alpha_grid), 1))
att_by_alpha <- matrix(NA, R_reps_path, length(alpha_grid))
catt_by_alpha <- array(NA, c(R_reps_path, length(alpha_grid), K))
att_by_m <- matrix(NA, R_reps_path, K); catt_by_m <- array(NA, c(R_reps_path, K, K))
for (r in seq_len(R_reps_path)) {
  d <- draw_rep(des_full, cohort_full, tau_true); c_vec <- d$c; yty <- d$yty
  tf <- as.vector(des_full$Ginv %*% c_vec)
  l0 <- l0_run(des_full, c_vec, yty, tf); msel_l0[r] <- l0$m_sel
  for (m in 1:K) { th <- l0$fits[[m]]$tau
    att_by_m[r, m] <- sum(w_att * th); catt_by_m[r, m, ] <- th }
  for (a in seq_along(alpha_grid)) {
    dp <- dp_bayes_rand(tf, c_vec, yty, des_full, alpha = alpha_grid[a])
    att_by_alpha[r, a] <- dp$att; catt_by_alpha[r, a, ] <- dp$tau
  }
}
## aggregate l0 path
l0_path <- data.frame(m = 1:K,
  att_bias = colMeans(att_by_m) - att_true,
  att_var  = apply(att_by_m, 2, var),
  catt_var = sapply(1:K, function(m) mean(apply(catt_by_m[, m, ], 2, var))))
## aggregate DP path
dp_path <- data.frame(alpha = alpha_grid,
  att_bias = colMeans(att_by_alpha) - att_true,
  att_var  = apply(att_by_alpha, 2, var),
  catt_var = sapply(seq_along(alpha_grid),
                    function(a) mean(apply(catt_by_alpha[, a, ], 2, var))))
write.csv(l0_path, file.path(out_dir, "l0_path.csv"), row.names = FALSE)
write.csv(dp_path, file.path(out_dir, "dp_path.csv"), row.names = FALSE)

png(file.path("figures", "solution_paths.png"), 2000, 900, res = 200)
par(mfrow = c(1, 2))
## l0 path vs m (reverse axis: more penalty -> fewer groups on the right)
plot(l0_path$m, abs(l0_path$att_bias), type = "b", pch = 19, col = "firebrick",
     xlab = "l0 groups m (small m = strong penalty)", ylab = "value",
     main = "l0 solution path (ATT bias & variance, CATT var)",
     ylim = range(0, abs(l0_path$att_bias), l0_path$att_var, l0_path$catt_var))
lines(l0_path$m, l0_path$att_var, type = "b", pch = 17, col = "steelblue")
lines(l0_path$m, l0_path$catt_var, type = "b", pch = 15, col = "darkgreen")
abline(v = mean(msel_l0), lty = 2)
legend("top", c("|ATT bias|","ATT var","CATT var","mean BIC m"),
       col = c("firebrick","steelblue","darkgreen","black"),
       pch = c(19,17,15,NA), lty = c(1,1,1,2), bty = "n")
## DP path vs alpha
plot(dp_path$alpha, abs(dp_path$att_bias), type = "b", pch = 19, col = "firebrick",
     xlab = expression(alpha~"(small="*strong~pooling*")"), ylab = "value",
     main = "DP solution path", log = "x",
     ylim = range(0, abs(dp_path$att_bias), dp_path$att_var, dp_path$catt_var))
lines(dp_path$alpha, dp_path$att_var, type = "b", pch = 17, col = "steelblue")
lines(dp_path$alpha, dp_path$catt_var, type = "b", pch = 15, col = "darkgreen")
legend("top", c("|ATT bias|","ATT var","CATT var"),
       col = c("firebrick","steelblue","darkgreen"), pch = c(19,17,15), bty = "n")
dev.off()
cat("PART 3 done -> l0_path.csv, dp_path.csv, solution_paths.png\n")

cat("\nAll parts complete.\n")
