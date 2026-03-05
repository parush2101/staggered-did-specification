# =============================================================================
# run_simulation.R — Main simulation script
# =============================================================================
#
# Demonstrates ℓ₀-penalized CATT estimation in a staggered DiD setting.
#
# Steps:
#   1. Generate panel data with a known CATT structure (2 true groups)
#   2. Run three estimators: Flexible, Pooled, L0-Greedy (multiple L values)
#   3. Report a summary table (# groups, RSS, RMSE)
#   4. Save plots to simulations/
# =============================================================================

# ── dependencies --------------------------------------------------------------
library(ggplot2)
source("R/dgp.R")
source("R/estimators.R")
source("R/plots.R")

set.seed(42)

# =============================================================================
# 1. Generate data
# =============================================================================
cat("=== Generating Panel Data ===\n")
panel    <- generate_panel(N_per_cohort = 100, T_periods = 8,
                           sigma = 0.5, seed = 42)
panel    <- create_treatment_dummies(panel)
panel_dm <- prepare_within(panel)

cat(sprintf("  Units: %d | Time periods: %d | Observations: %d\n",
            length(unique(panel$unit)),
            length(unique(panel$time)),
            nrow(panel)))
cat("  Cohort sizes (units):\n")
print(table(panel$cohort[panel$time == 1]))

catts <- get_catt_list(max(panel$time))
cat(sprintf("\n  Total CATTs: %d\n", nrow(catts)))
cat("  True partition: 2 groups\n")
cat("    Group A (τ = 2): cohorts 3 and 5 post-treatment cells (10 CATTs)\n")
cat("    Group B (τ = 4): cohort 7 post-treatment cells          (2 CATTs)\n\n")

# =============================================================================
# 2. Fully Flexible
# =============================================================================
cat("=== Fully Flexible Estimator ===\n")
flex_res    <- estimate_flexible(panel_dm)
flex_ests   <- collect_estimates(flex_res, panel_dm, "Flexible")
flex_rmse   <- sqrt(mean((flex_ests$estimate - flex_ests$tau_true)^2))
cat(sprintf("  Groups: %d | RSS: %.4f | RMSE: %.4f\n\n",
            length(flex_res$partition), flex_res$rss, flex_rmse))

# =============================================================================
# 3. Fully Pooled
# =============================================================================
cat("=== Fully Pooled Estimator ===\n")
pool_res    <- estimate_pooled(panel_dm)
pool_ests   <- collect_estimates(pool_res, panel_dm, "Pooled")
pool_rmse   <- sqrt(mean((pool_ests$estimate - pool_ests$tau_true)^2))
cat(sprintf("  Groups: %d | RSS: %.4f | RMSE: %.4f\n\n",
            length(pool_res$partition), pool_res$rss, pool_rmse))

# =============================================================================
# 4. ℓ₀ Greedy Estimator — multiple values of L
# =============================================================================
# Intuition for L choices:
#   Small L   → accepts almost any merge → risks over-pooling
#   Large L   → accepts only very cheap merges → risks under-pooling
#   Correct L → recovers the true 2-group partition
L_values <- c(0.5, 2, 10, 50)

cat("=== L0 Greedy Estimator ===\n")
l0_results <- lapply(L_values, function(L) {
  cat(sprintf("-- L = %.1f --\n", L))
  greedy_l0(panel_dm, L = L, verbose = TRUE)
})
names(l0_results) <- paste0("L0 (L=", L_values, ")")

# =============================================================================
# 5. Summary table
# =============================================================================
l0_ests_list <- mapply(
  collect_estimates,
  result      = l0_results,
  method_label = names(l0_results),
  MoreArgs    = list(panel_dm = panel_dm),
  SIMPLIFY    = FALSE
)

make_row <- function(label, result, ests) {
  data.frame(
    Method   = label,
    N_Groups = length(result$partition),
    RSS      = round(result$rss, 4),
    RMSE     = round(sqrt(mean((ests$estimate - ests$tau_true)^2)), 4),
    stringsAsFactors = FALSE
  )
}

summary_df <- rbind(
  make_row("Flexible", flex_res, flex_ests),
  do.call(rbind, mapply(make_row,
                        label  = names(l0_results),
                        result = l0_results,
                        ests   = l0_ests_list,
                        SIMPLIFY = FALSE)),
  make_row("Pooled", pool_res, pool_ests)
)

cat("\n=== Summary Table ===\n")
print(summary_df, row.names = FALSE)

# =============================================================================
# 6. Plots
# =============================================================================

# 6a. CATT estimates: Flexible, L0 (L=2), Pooled
# L=2 is expected to recover the true 2-group partition in this DGP
target_l0 <- "L0 (L=2)"

plot_df <- rbind(flex_ests, l0_ests_list[[target_l0]], pool_ests)
p_catts  <- plot_catts(plot_df)

ggsave("simulations/catt_estimates.png", p_catts,
       width = 10, height = 4, dpi = 150)
cat("\nSaved: simulations/catt_estimates.png\n")

# 6b. Solution path: number of groups vs L
path_df <- data.frame(
  L        = L_values,
  n_groups = sapply(l0_results, `[[`, "n_groups"),
  rss      = sapply(l0_results, `[[`, "rss")
)

p_path <- plot_solution_path(path_df)
ggsave("simulations/solution_path.png", p_path,
       width = 7, height = 4, dpi = 150)
cat("Saved: simulations/solution_path.png\n")

cat("\nDone.\n")
