# =============================================================================
# monte_carlo.R — Monte Carlo evaluation of CATT estimators
# =============================================================================
#
# Runs S = 1000 independent replications of the DGP, applies all estimators,
# and summarises the sampling distributions of the estimated ATTs.
#
# Output:
#   - Console: mean and variance table by method and true group
#   - simulations/mc_summary.csv       — group-level summary table
#   - simulations/mc_catt_estimates.csv — CATT-level mean, variance, bias
#   - simulations/mc_sampling_dist.png  — sampling distributions by group
#   - simulations/mc_catt_variance.png  — variance of each CATT by method
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

# =============================================================================
# 4. CATT-level summary: mean and variance for each (method, cohort, time) cell
# =============================================================================
catt_mean <- aggregate(estimate ~ method + cohort + time + label + tau_true,
                       data = mc_df, FUN = mean)
catt_var  <- aggregate(estimate ~ method + cohort + time + label + tau_true,
                       data = mc_df, FUN = var)

catt_df           <- catt_mean
names(catt_df)[names(catt_df) == "estimate"] <- "mean_est"
catt_df$var_est   <- catt_var$estimate
catt_df$bias      <- round(catt_df$mean_est - catt_df$tau_true, 6)
catt_df$mean_est  <- round(catt_df$mean_est, 6)
catt_df$var_est   <- round(catt_df$var_est,  8)
catt_df$method    <- factor(catt_df$method, levels = method_order)
catt_df           <- catt_df[order(catt_df$method, catt_df$cohort, catt_df$time), ]
rownames(catt_df) <- NULL

write.csv(catt_df, "simulations/mc_catt_estimates.csv", row.names = FALSE)
cat("Saved: simulations/mc_catt_estimates.csv\n")

# =============================================================================
# 5. Variance plot: one line per method, x = CATT, y = variance of estimate
# =============================================================================

# Order CATTs by cohort then time and build readable x-axis labels
catt_order  <- with(catt_df[catt_df$method == "Flexible", ],
                    label[order(cohort, time)])
catt_labels <- gsub("D_", "C", catt_order)          # D_3_5 -> C3_5
catt_labels <- gsub("_", ",t=", catt_labels)         # C3_5  -> C3,t=5

catt_df$label_f <- factor(catt_df$label, levels = catt_order)

# Cohort bands for background shading
band_df <- data.frame(
  cohort = c(3, 5, 7),
  xmin   = c(0.5, 6.5, 10.5),
  xmax   = c(6.5, 10.5, 12.5),
  fill   = c("grey90", "grey80", "grey70")
)

p_var <- ggplot(catt_df, aes(x = label_f, y = var_est,
                              color = method, group = method)) +
  annotate("rect", xmin = band_df$xmin, xmax = band_df$xmax,
           ymin = -Inf, ymax = Inf, fill = band_df$fill, alpha = 0.3) +
  annotate("text",
           x     = c(3.5, 8.5, 11.5),
           y     = Inf, vjust = 1.5, size = 3.2, color = "grey30",
           label = c("Cohort 3 (tau=2)", "Cohort 5 (tau=2)", "Cohort 7 (tau=4)")) +
  geom_line(linewidth = 0.85) +
  geom_point(size = 2.2) +
  scale_x_discrete(labels = catt_labels) +
  scale_color_manual(values = c(
    "Flexible"    = "#E41A1C",
    "L0 (L=0.5)"  = "#C6DBEF",
    "L0 (L=2)"    = "#6BAED6",
    "L0 (L=10)"   = "#2171B5",
    "L0 (L=50)"   = "#08306B",
    "Pooled"      = "#4DAF4A"
  )) +
  labs(
    title    = expression(paste("Variance of CATT Estimates Across ", italic(S), " = 1000 Replications")),
    subtitle = "Lower variance = more efficient estimator",
    x        = "CATT  (cohort, calendar time)",
    y        = "Variance of Estimated CATT",
    color    = "Method"
  ) +
  theme_bw(base_size = 12) +
  theme(
    axis.text.x      = element_text(angle = 45, hjust = 1, size = 9),
    legend.position  = "bottom",
    panel.grid.minor = element_blank()
  )

ggsave("simulations/mc_catt_variance.png", p_var,
       width = 11, height = 5, dpi = 150)
cat("Saved: simulations/mc_catt_variance.png\n")

cat("\nDone.\n")
