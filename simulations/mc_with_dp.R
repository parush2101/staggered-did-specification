# =============================================================================
# mc_with_dp.R — Monte Carlo with all four estimators incl. DP Bayes
# =============================================================================
#
# S = 500 replications.  For each draw:
#   - Flexible, L0 (L=0.5), Pooled: point estimates as before
#   - DP Bayes: 100 Gibbs iterations; posterior mean per CATT
#
# Outputs:
#   simulations/mc_dp_summary.csv   — group-level table for the paper
#   simulations/mc_dp_iter_means.csv — per-replication group means (all methods)
# =============================================================================

library(ggplot2)
source("R/dgp.R")
source("R/estimators.R")

S      <- 500
n_iter <- 100          # Gibbs iterations per replication
alpha    <- 1
sigma2_0 <- 9
mu_0     <- 0

set.seed(123)
seeds <- sample.int(1e6, S)

cat(sprintf("Running %d MC replications  (DP Gibbs: %d iters each)...\n",
            S, n_iter))
cat("Progress: ")

mc_list <- vector("list", S)

for (s in seq_len(S)) {

  if (s %% 50 == 0) cat(sprintf("%d ", s))

  # ---- Data ----
  panel    <- generate_panel(N_per_cohort = 100, T_periods = 8,
                             sigma = 0.5, seed = seeds[s])
  panel    <- create_treatment_dummies(panel)
  panel_dm <- prepare_within(panel)

  catts <- get_catt_list(max(panel_dm$time))
  K     <- nrow(catts)
  NT    <- nrow(panel_dm)
  y_dm  <- panel_dm$y_dm

  # ---- Flexible ----
  flex_res  <- estimate_flexible(panel_dm)
  flex_ests <- collect_estimates(flex_res, panel_dm, "Flexible")

  # ---- Pooled ----
  pool_res  <- estimate_pooled(panel_dm)
  pool_ests <- collect_estimates(pool_res, panel_dm, "Pooled")

  # ---- L0 (L=0.5) ----
  l0_res  <- greedy_l0(panel_dm, L = 0.5, verbose = FALSE)
  l0_ests <- collect_estimates(l0_res, panel_dm, "L0 (L=0.5)")

  # ---- DP Bayes ----
  N_units   <- length(unique(panel_dm$unit))
  T_periods <- length(unique(panel_dm$time))
  df        <- NT - (N_units + T_periods - 1) - K
  sigma2    <- flex_res$rss / df
  kappa     <- sigma2_0 / sigma2

  # Log marginal likelihood (closed form after integrating out phi)
  log_marg_lik <- function(partition) {
    m <- length(partition)
    X <- matrix(0, nrow = NT, ncol = m)
    for (l in seq_len(m)) {
      dm_cols <- paste0(partition[[l]], "_dm")
      X[, l]  <- rowSums(panel_dm[, dm_cols, drop = FALSE])
    }
    A   <- diag(m) + kappa * crossprod(X)
    Xty <- crossprod(X, y_dm)
    -0.5 * as.numeric(determinant(A, logarithm = TRUE)$modulus) +
      kappa / (2 * sigma2) * as.numeric(t(Xty) %*% solve(A) %*% Xty)
  }

  # Collapsed Gibbs sampler
  z           <- seq_len(K)   # initialise: each CATT its own singleton
  theta_draws <- matrix(NA_real_, nrow = n_iter, ncol = K)

  for (iter in seq_len(n_iter)) {

    # Update each CATT's cluster assignment
    for (k in seq_len(K)) {
      z_minus_k       <- z[-k]
      unique_clusters <- sort(unique(z_minus_k))
      n_minus_k       <- tabulate(z_minus_k)
      c_new           <- max(z) + 1L
      candidates      <- c(unique_clusters, c_new)
      log_probs       <- numeric(length(candidates))

      for (idx in seq_along(candidates)) {
        c_cand     <- candidates[idx]
        z_trial    <- z; z_trial[k] <- c_cand
        lbl_map    <- sort(unique(z_trial))
        z_rel      <- match(z_trial, lbl_map)
        partition  <- lapply(seq_along(lbl_map),
                             function(l) catts$label[z_rel == l])
        crp_wt     <- if (c_cand == c_new) alpha else n_minus_k[c_cand]
        log_probs[idx] <- log(crp_wt) + log_marg_lik(partition)
      }

      log_probs <- log_probs - max(log_probs)
      probs     <- exp(log_probs) / sum(exp(log_probs))
      z[k]      <- candidates[sample(length(candidates), 1L, prob = probs)]
      z         <- match(z, sort(unique(z)))
    }

    # Draw group effects phi from conjugate posterior
    m_curr    <- max(z)
    partition <- lapply(seq_len(m_curr), function(l) catts$label[z == l])
    X <- matrix(0, nrow = NT, ncol = m_curr)
    for (l in seq_len(m_curr)) {
      dm_cols <- paste0(partition[[l]], "_dm")
      X[, l]  <- rowSums(panel_dm[, dm_cols, drop = FALSE])
    }
    Sigma_star <- solve(crossprod(X) / sigma2 + diag(m_curr) / sigma2_0)
    mu_star    <- Sigma_star %*%
      (crossprod(X, y_dm) / sigma2 + mu_0 * rep(1, m_curr) / sigma2_0)
    L_chol     <- tryCatch(
      chol(Sigma_star),
      error = function(e) chol(Sigma_star + 1e-8 * diag(m_curr))
    )
    phi <- as.numeric(mu_star) + as.numeric(L_chol %*% rnorm(m_curr))
    theta_draws[iter, ] <- phi[z]
  }

  # Posterior mean for each CATT across all iterations
  dp_post_means <- colMeans(theta_draws)

  dp_ests <- data.frame(
    cohort   = catts$cohort,
    time     = catts$time,
    label    = catts$label,
    tau_true = catts$tau_true,
    estimate = dp_post_means,
    group    = NA_integer_,
    method   = "DP Bayes",
    stringsAsFactors = FALSE
  )

  iter_df           <- rbind(flex_ests, pool_ests, l0_ests, dp_ests)
  iter_df$iteration <- s
  mc_list[[s]]      <- iter_df
}

cat("\nDone.\n\n")

mc_df <- do.call(rbind, mc_list)

# True group labels
mc_df$true_group <- ifelse(mc_df$tau_true == 2,
                            "Group A (tau=2)", "Group B (tau=4)")

# Per-replication group means
iter_means <- aggregate(
  estimate ~ method + true_group + tau_true + iteration,
  data = mc_df,
  FUN  = mean
)

write.csv(iter_means, "simulations/mc_dp_iter_means.csv", row.names = FALSE)

# Summary table
method_order <- c("Flexible", "L0 (L=0.5)", "DP Bayes", "Pooled")

summary_df <- do.call(rbind, lapply(
  split(iter_means,
        list(iter_means$method, iter_means$true_group),
        drop = TRUE),
  function(sub) {
    data.frame(
      Method     = sub$method[1],
      True_Group = sub$true_group[1],
      True_tau   = sub$tau_true[1],
      Mean_Est   = round(mean(sub$estimate), 4),
      Var_Est    = round(var(sub$estimate),  6),
      Bias       = round(mean(sub$estimate) - sub$tau_true[1], 4),
      stringsAsFactors = FALSE
    )
  }
))

summary_df$Method <- factor(summary_df$Method, levels = method_order)
summary_df <- summary_df[order(summary_df$Method, summary_df$True_Group), ]
rownames(summary_df) <- NULL

cat("=== Monte Carlo Summary Table (S = 500) ===\n")
cat("Mean_Est = mean of group-level ATT estimates across replications\n")
cat("Var_Est  = variance of group-level ATT estimates across replications\n")
cat("Bias     = Mean_Est - true tau\n\n")
print(summary_df, row.names = FALSE)

write.csv(summary_df, "simulations/mc_dp_summary.csv", row.names = FALSE)
cat("\nSaved: simulations/mc_dp_summary.csv\n")
cat("Saved: simulations/mc_dp_iter_means.csv\n")
cat("\nDone.\n")
