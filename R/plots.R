# =============================================================================
# plots.R — Visualisation for CATT simulation results
# =============================================================================

library(ggplot2)

# -----------------------------------------------------------------------------
# plot_catts()
#
# Line plot of CATT estimates over calendar time, faceted by cohort.
# The true effect is shown as a black dashed horizontal line.
#
# Arguments:
#   df — data frame produced by rbind-ing collect_estimates() outputs
# -----------------------------------------------------------------------------
plot_catts <- function(df) {

  # Factor levels control legend order
  df$method <- factor(df$method,
                      levels = c("Flexible", "Pooled",
                                 grep("^L0", unique(df$method), value = TRUE)))

  method_colors <- c(
    "Flexible" = "#E41A1C",
    "Pooled"   = "#4DAF4A"
  )
  # Assign blues to L0 variants
  l0_methods <- grep("^L0", levels(df$method), value = TRUE)
  blues       <- colorRampPalette(c("#9ECAE1", "#08519C"))(length(l0_methods))
  names(blues) <- l0_methods
  all_colors   <- c(method_colors, blues)

  ggplot(df, aes(x = time, y = estimate,
                 color = method, group = method)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 2.2) +
    # True effect as horizontal dashed reference
    geom_segment(
      aes(x = min(time) - 0.3, xend = max(time) + 0.3,
          y = tau_true,         yend = tau_true),
      color     = "black",
      linetype  = "dashed",
      linewidth = 0.8,
      inherit.aes = FALSE,
      data = df[!duplicated(df[, c("cohort", "tau_true")]), ]
    ) +
    facet_wrap(~ paste("Cohort", cohort), ncol = 3) +
    scale_color_manual(values = all_colors) +
    scale_x_continuous(breaks = scales::pretty_breaks()) +
    labs(
      title    = expression(paste(ell[0], "-Penalized CATT Estimation")),
      subtitle = "Black dashed line = true effect",
      x        = "Calendar Time",
      y        = "Estimated CATT",
      color    = "Method"
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position  = "bottom",
      strip.background = element_rect(fill = "grey92"),
      panel.grid.minor = element_blank()
    )
}


# -----------------------------------------------------------------------------
# plot_solution_path()
#
# Shows how the number of groups and RSS change as L increases.
#
# Arguments:
#   summary_df — data frame with columns: L, n_groups, rss
# -----------------------------------------------------------------------------
plot_solution_path <- function(summary_df) {

  p1 <- ggplot(summary_df, aes(x = log10(L), y = n_groups)) +
    geom_line(linewidth = 0.9, color = "#2171B5") +
    geom_point(size = 3, color = "#2171B5") +
    geom_hline(yintercept = 2, linetype = "dashed", color = "grey40") +
    annotate("text", x = min(log10(summary_df$L)) + 0.1,
             y = 2.15, label = "True # groups = 2",
             hjust = 0, size = 3.5, color = "grey40") +
    scale_x_continuous(
      breaks = log10(summary_df$L),
      labels = summary_df$L
    ) +
    labs(
      title = "Solution Path: Number of Groups vs. Penalty L",
      x     = "Penalty L (log scale)",
      y     = "Number of Groups"
    ) +
    theme_bw(base_size = 12) +
    theme(panel.grid.minor = element_blank())

  p1
}
