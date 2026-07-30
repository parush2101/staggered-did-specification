# =============================================================================
# dp_posterior.R — Bayesian CATT estimation via Dirichlet Process mixture
# =============================================================================
#
# Runs a collapsed Gibbs sampler (Neal 2000, Algorithm 3) on the original DGP.
# The DP prior clusters CATTs that share the same treatment effect.
#
# Steps:
#   1. Generate data (seed = 42, same as run_simulation.R)
#   2. Estimate hyperparameters from the flexible model
#   3. Run 100 Gibbs iterations
#   4. Plot posterior densities of CATT effects, faceted by cohort
#
# Outputs:
#   simulations/dp_posterior.png
# =============================================================================

library(ggplot2)
source("R/dgp.R")
source("R/estimators.R")

# =============================================================================
# 1. Data
# =============================================================================
set.seed(42)
panel    <- generate_panel(N_per_cohort = 100, T_periods = 8, sigma = 0.5, seed = 42)
panel    <- create_treatment_dummies(panel)
panel_dm <- prepare_within(panel)

catts <- get_catt_list(max(panel_dm$time))
K     <- nrow(catts)
NT    <- nrow(panel_dm)
y_dm  <- panel_dm$y_dm

cat(sprintf("Data: NT = %d, K = %d CATTs\n", NT, K))

# =============================================================================
# 2. Hyperparameters
# =============================================================================
flex_res <- estimate_flexible(panel_dm)

# Noise variance: residual from flexible model
N_units   <- length(unique(panel_dm$unit))
T_periods <- length(unique(panel_dm$time))
df        <- NT - (N_units + T_periods - 1) - K
sigma2    <- flex_res$rss / df

# Base measure: N(mu_0, sigma2_0)
# Centre at zero; prior SD covers the range of treatment effects
mu_0     <- 0
sigma2_0 <- 9          # prior SD = 3 (covers effects 0-6 within 2 SDs)
kappa    <- sigma2_0 / sigma2

# DP concentration: alpha = 1 gives ~log(K) ~ 2.5 expected clusters
alpha <- 1

cat(sprintf("sigma2 = %.4f  |  sigma2_0 = %.1f  |  alpha = %.1f\n\n",
            sigma2, sigma2_0, alpha))

# =============================================================================
# 3. Log marginal likelihood (integrating out group effects phi)
#
# p(y_dm | P) = N(0, sigma2*I + sigma2_0 * X_tilde * X_tilde')
#
# Terms varying across partitions (constant terms dropped):
#   log p(y_dm | P) ∝  -1/2 log|I_m + kappa * X'X|
#                      + kappa/(2*sigma2) * (X'y)' (I_m + kappa*X'X)^{-1} (X'y)
# =============================================================================
log_marg_lik <- function(partition) {
  m <- length(partition)

  # Group regressors: X_tilde_l = sum of within-transformed dummies in group l
  X <- matrix(0, nrow = NT, ncol = m)
  for (l in seq_len(m)) {
    dm_cols <- paste0(partition[[l]], "_dm")
    X[, l]  <- rowSums(panel_dm[, dm_cols, drop = FALSE])
  }

  XtX  <- crossprod(X)                       # m x m
  A    <- diag(m) + kappa * XtX              # I_m + kappa * X'X
  Xty  <- crossprod(X, y_dm)                 # m x 1

  A_inv   <- solve(A)
  log_det <- as.numeric(determinant(A, logarithm = TRUE)$modulus)
  quad    <- as.numeric(t(Xty) %*% A_inv %*% Xty)

  -0.5 * log_det + kappa / (2 * sigma2) * quad
}

# =============================================================================
# 4. Collapsed Gibbs sampler
# =============================================================================
n_iter <- 100
z      <- seq_len(K)           # initialise: each CATT in its own singleton cluster

theta_draws <- matrix(NA_real_, nrow = n_iter, ncol = K)
colnames(theta_draws) <- catts$label
m_draws     <- integer(n_iter)

cat(sprintf("Running %d Gibbs iterations...\n", n_iter))

for (iter in seq_len(n_iter)) {

  # ---- update each CATT's cluster assignment --------------------------------
  for (k in seq_len(K)) {

    z_minus_k       <- z[-k]
    unique_clusters <- sort(unique(z_minus_k))
    n_minus_k       <- tabulate(z_minus_k)   # counts excl. k

    # Candidates: all existing non-empty clusters + one new cluster
    c_new      <- max(z) + 1L
    candidates <- c(unique_clusters, c_new)
    log_probs  <- numeric(length(candidates))

    for (idx in seq_along(candidates)) {
      c_cand    <- candidates[idx]
      z_trial   <- z
      z_trial[k] <- c_cand

      # Re-label clusters consecutively
      lbl_map   <- sort(unique(z_trial))
      z_rel     <- match(z_trial, lbl_map)
      partition <- lapply(seq_along(lbl_map),
                          function(l) catts$label[z_rel == l])

      # CRP weight
      crp_wt <- if (c_cand == c_new) alpha else n_minus_k[c_cand]

      log_probs[idx] <- log(crp_wt) + log_marg_lik(partition)
    }

    # Normalise and sample
    log_probs <- log_probs - max(log_probs)
    probs     <- exp(log_probs); probs <- probs / sum(probs)
    z[k]      <- candidates[sample(length(candidates), 1L, prob = probs)]

    # Keep labels consecutive
    lbl  <- sort(unique(z))
    z    <- match(z, lbl)
  }

  # ---- draw group effects phi from their posterior --------------------------
  m_curr    <- max(z)
  partition <- lapply(seq_len(m_curr), function(l) catts$label[z == l])

  X <- matrix(0, nrow = NT, ncol = m_curr)
  for (l in seq_len(m_curr)) {
    dm_cols <- paste0(partition[[l]], "_dm")
    X[, l]  <- rowSums(panel_dm[, dm_cols, drop = FALSE])
  }

  Sigma_star <- solve(crossprod(X) / sigma2 + diag(m_curr) / sigma2_0)
  mu_star    <- Sigma_star %*% (crossprod(X, y_dm) / sigma2 +
                                 mu_0 * rep(1, m_curr) / sigma2_0)

  # Cholesky draw
  L_chol <- tryCatch(
    chol(Sigma_star),
    error = function(e) chol(Sigma_star + 1e-8 * diag(m_curr))
  )
  phi <- as.numeric(mu_star) + as.numeric(L_chol %*% rnorm(m_curr))

  theta_draws[iter, ] <- phi[z]
  m_draws[iter]       <- m_curr

  if (iter %% 20 == 0)
    cat(sprintf("  Iter %3d | groups: %d\n", iter, m_curr))
}

cat(sprintf(
  "\nPosterior # groups: mean = %.1f | min = %d | max = %d\n\n",
  mean(m_draws), min(m_draws), max(m_draws)
))

# =============================================================================
# 5. Posterior co-clustering probabilities
# =============================================================================
co_clust <- matrix(0, K, K, dimnames = list(catts$label, catts$label))
for (iter in seq_len(n_iter)) {
  z_iter <- match(theta_draws[iter, ], unique(theta_draws[iter, ]))
  for (j in seq_len(K))
    for (k in seq_len(K))
      if (z_iter[j] == z_iter[k]) co_clust[j, k] <- co_clust[j, k] + 1
}
co_clust <- co_clust / n_iter

cat("=== Posterior co-clustering probability matrix (rows/cols = CATTs) ===\n")
cat("(Values close to 1 = almost always in the same group)\n\n")
grp_A <- catts$label[catts$tau_true == 2]
grp_B <- catts$label[catts$tau_true == 4]
cat(sprintf("Mean P(same group | both in Group A, tau=2): %.3f\n",
            mean(co_clust[grp_A, grp_A][lower.tri(co_clust[grp_A, grp_A])])))
cat(sprintf("Mean P(same group | both in Group B, tau=4): %.3f\n",
            mean(co_clust[grp_B, grp_B][lower.tri(co_clust[grp_B, grp_B])])))
cat(sprintf("Mean P(same group | one in A, one in B):     %.3f\n\n",
            mean(co_clust[grp_A, grp_B])))

# =============================================================================
# 6. Plot posterior densities
# =============================================================================
draws_long <- do.call(rbind, lapply(seq_len(K), function(k) {
  data.frame(
    label    = catts$label[k],
    cohort   = catts$cohort[k],
    time     = catts$time[k],
    tau_true = catts$tau_true[k],
    theta    = theta_draws[, k]
  )
}))

draws_long$time_f <- factor(draws_long$time)

# True effect reference (one row per cohort)
ref_df <- unique(draws_long[, c("cohort", "tau_true")])

p <- ggplot(draws_long,
            aes(x = theta, color = time_f, fill = time_f, group = label)) +
  geom_density(alpha = 0.15, linewidth = 0.7) +
  geom_vline(data = ref_df,
             aes(xintercept = tau_true),
             linetype = "dashed", color = "black", linewidth = 0.8,
             inherit.aes = FALSE) +
  facet_wrap(~ paste("Cohort", cohort), ncol = 3, scales = "free_y") +
  scale_color_brewer(palette = "Dark2") +
  scale_fill_brewer(palette  = "Dark2") +
  labs(
    title    = "Posterior Distributions of CATT Effects — DP Mixture Gibbs Sampler",
    subtitle = sprintf(
      "%d iterations  |  K = %d CATTs  |  \u03b1 = %.1f  |  Black dashed = true effect",
      n_iter, K, alpha),
    x        = "CATT Effect",
    y        = "Density",
    color    = "Calendar Time",
    fill     = "Calendar Time"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "bottom",
    strip.background = element_rect(fill = "grey92"),
    panel.grid.minor = element_blank()
  )

ggsave("simulations/dp_posterior.png", p, width = 11, height = 5, dpi = 150)
cat("Saved: simulations/dp_posterior.png\n")
