# =============================================================================
# monte_carlo.R — Monte Carlo evaluation of CATT estimators
# =============================================================================
#
# Runs S = 1000 independent replications of the DGP, applies all estimators,
# and summarises the sampling distributions of the estimated ATTs.
#
# Output:
#   - Console: mean and variance table by method and true group
#   - simulations/mc_summary.csv  — full summary table
#   - simulations/mc_sampling_dist.png — sampling distributions plot
# =============================================================================

library(ggplot2)
source("R/dgp.R")
source("R/estimators.R")

S        <- 1000
L_values <- c(0.5, 2, 10, 50)
set.seed(123)
seeds <- sample.int(1e6, S)

cat(sprintf("Running %d Monte Carlo replications...\n", S))
cat("Progress: ")

# =============================================================================
# 1. Run replications
# =============================================================================
mc_list <- vector("list", S)

for (s in seq_len(S)) {

  if (s %% 100 == 0) cat(sprintf("%d ", s))

  panel    <- generate_panel(N_per_cohort = 100, T_periods = 8,
                             sigma = 0.5, seed = seeds[s])
  panel    <- create_treatment_dummies(panel)
  panel_dm <- prepare_within(panel)

  # Flexible
  flex_res  <- estimate_flexible(panel_dm)
  flex_ests <- collect_estimates(flex_res, panel_dm, "Flexible")

  # Pooled
  pool_res  <- estimate_pooled(panel_dm)
  pool_ests <- collect_estimates(pool_res, panel_dm, "Pooled")

  # L0 Greedy for each L
  l0_ests <- lapply(L_values, function(L) {
    res <- greedy_l0(panel_dm, L = L, verbose = FALSE)
    collect_estimates(res, panel_dm, paste0("L0 (L=", L, ")"))
  })

  iter_df <- rbind(flex_ests, pool_ests, do.call(rbind, l0_ests))
  iter_df$iteration <- s
  mc_list[[s]] <- iter_df
}

cat("\nDone.\n\n")

mc_df <- do.call(rbind, mc_list)

# =============================================================================
# 2. Summary table: mean and variance by method × true group
# =============================================================================
# True groups: Group A (tau = 2), Group B (tau = 4)
mc_df$true_group <- ifelse(mc_df$tau_true == 2, "Group A (τ=2)", "Group B (τ=4)")

# Collapse to method × true_group × iteration first (average within group per draw)
iter_means <- aggregate(
  estimate ~ method + true_group + tau_true + iteration,
  data = mc_df,
  FUN  = mean
)

# Then compute mean and variance across iterations
summary_df <- do.call(rbind, lapply(
  split(iter_means, list(iter_means$method, iter_means$true_group), drop = TRUE),
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

# Order rows: Flexible → L0 variants (ascending L) → Pooled
method_order <- c("Flexible",
                  paste0("L0 (L=", L_values, ")"),
                  "Pooled")
summary_df$Method <- factor(summary_df$Method, levels = method_order)
summary_df <- summary_df[order(summary_df$Method, summary_df$True_Group), ]
rownames(summary_df) <- NULL

cat("=== Monte Carlo Summary Table (S = 1000) ===\n")
cat("Columns: Mean_Est = mean of estimated ATT across replications\n")
cat("         Var_Est  = variance of estimated ATT across replications\n")
cat("         Bias     = Mean_Est - true tau\n\n")
print(summary_df, row.names = FALSE)

# Save to CSV
write.csv(summary_df, "simulations/mc_summary.csv", row.names = FALSE)
cat("\nSaved: simulations/mc_summary.csv\n")

# =============================================================================
# 3. Sampling distribution plot
# =============================================================================
# Show density of iteration-level group means for each method
plot_df <- iter_means
plot_df$method <- factor(plot_df$method, levels = method_order)

p <- ggplot(plot_df, aes(x = estimate, fill = method, color = method)) +
  geom_density(alpha = 0.25, linewidth = 0.6) +
  geom_vline(aes(xintercept = tau_true), linetype = "dashed",
             color = "black", linewidth = 0.7) +
  facet_wrap(~ true_group, scales = "free") +
  labs(
    title    = expression(paste("Sampling Distributions of Estimated ATT (", italic(S), " = 1000)")),
    subtitle = "Black dashed line = true effect",
    x        = "Estimated ATT",
    y        = "Density",
    fill     = "Method",
    color    = "Method"
  ) +
  theme_bw(base_size = 12) +
  theme(
    legend.position  = "bottom",
    strip.background = element_rect(fill = "grey92"),
    panel.grid.minor = element_blank()
  )

ggsave("simulations/mc_sampling_dist.png", p,
       width = 10, height = 5, dpi = 150)
cat("Saved: simulations/mc_sampling_dist.png\n")

cat("\nDone.\n")
