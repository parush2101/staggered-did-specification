# =============================================================================
# dgp.R — Data Generating Process for Staggered DiD Simulation
# =============================================================================
#
# Panel structure:
#   - N = 400 units (100 per cohort × 4 cohorts)
#   - T = 8 time periods
#   - Three treated cohorts (enter treatment at t = 3, 5, 7)
#   - One never-treated cohort
#
# True CATT structure (the target our estimator should recover):
#   - Cohort 3 → τ = 2  for t = 3, …, 8   (6 CATTs)
#   - Cohort 5 → τ = 2  for t = 5, …, 8   (4 CATTs)   ← same as cohort 3
#   - Cohort 7 → τ = 4  for t = 7, 8       (2 CATTs)   ← different
#
# True partition: 2 groups
#   Group A (10 CATTs): all cohort-3 and cohort-5 post-treatment cells, τ = 2
#   Group B ( 2 CATTs): all cohort-7 post-treatment cells,              τ = 4
#
# Outcome model:
#   y_it = α_i + λ_t + τ_gt · D_git + ε_it
#   α_i  ~ N(0, 1)          (unit fixed effect)
#   λ_t  = 0.1 · (t - 1)   (linear time trend)
#   ε_it ~ N(0, σ²)
# =============================================================================

generate_panel <- function(N_per_cohort = 100,
                           T_periods    = 8,
                           sigma        = 0.5,
                           seed         = 42) {
  set.seed(seed)

  # --- cohort definitions ---------------------------------------------------
  cohort_entry <- c(3, 5, 7, Inf)          # Inf = never treated
  true_tau     <- c("3" = 2, "5" = 2, "7" = 4)

  # --- unit-level data ------------------------------------------------------
  N      <- N_per_cohort * length(cohort_entry)
  units  <- data.frame(
    unit    = seq_len(N),
    cohort  = rep(cohort_entry, each = N_per_cohort),
    alpha_i = rnorm(N, mean = 0, sd = 1)
  )

  # --- expand to balanced panel ---------------------------------------------
  panel <- merge(units, data.frame(time = seq_len(T_periods)), by = NULL)
  panel <- panel[order(panel$unit, panel$time), ]
  rownames(panel) <- NULL

  # --- time fixed effect (deterministic linear trend) -----------------------
  panel$lambda_t <- 0.1 * (panel$time - 1)

  # --- treatment indicator --------------------------------------------------
  panel$w_it <- as.integer(
    !is.infinite(panel$cohort) & panel$time >= panel$cohort
  )

  # --- true treatment effect ------------------------------------------------
  panel$tau_true <- ifelse(
    panel$w_it == 1,
    true_tau[as.character(panel$cohort)],
    0
  )

  # --- outcome --------------------------------------------------------------
  panel$y <- panel$alpha_i + panel$lambda_t + panel$tau_true +
             rnorm(nrow(panel), mean = 0, sd = sigma)

  return(panel)
}


# -----------------------------------------------------------------------------
# get_catt_list()
#
# Returns a data frame of all cohort-time cells that correspond to a CATT,
# i.e., all (g, t) pairs where t >= g and cohort g is treated.
#
# Columns: cohort, time, label ("D_g_t"), tau_true
# -----------------------------------------------------------------------------
get_catt_list <- function(T_periods = 8) {
  true_tau <- c("3" = 2, "5" = 2, "7" = 4)
  treated_cohorts <- c(3, 5, 7)

  do.call(rbind, lapply(treated_cohorts, function(g) {
    ts <- g:T_periods
    data.frame(
      cohort   = g,
      time     = ts,
      label    = paste0("D_", g, "_", ts),
      tau_true = unname(true_tau[as.character(g)]),
      stringsAsFactors = FALSE
    )
  }))
}


# -----------------------------------------------------------------------------
# create_treatment_dummies()
#
# Adds a binary column D_g_t for every CATT (g, t) in the panel.
# D_g_t = 1 iff unit belongs to cohort g AND current time period == t.
# -----------------------------------------------------------------------------
create_treatment_dummies <- function(panel) {
  catts <- get_catt_list(max(panel$time))

  for (i in seq_len(nrow(catts))) {
    g   <- catts$cohort[i]
    t   <- catts$time[i]
    lbl <- catts$label[i]
    panel[[lbl]] <- as.integer(panel$cohort == g & panel$time == t)
  }
  return(panel)
}
