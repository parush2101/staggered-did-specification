# ============================================================
# Monte Carlo Simulation for Staggered DiD Specification
#
# Design grid:
#   - T = 10
#   - 4 cohorts: never-treated + treated at t = 3, 5, 7
#   - Cohort sizes: small only - uniform_small (all 10) and unequal_small (5/10/13/15).
#                   Homogeneous truth is run with uniform_small only (cohort imbalance
#                   carries no extra information when CATTs are constant).
#   - CATT structures: homogeneous, 6 partitions (mixed-gap), fully heterogeneous
#   - 100 Monte Carlo replications per scenario
#
# Estimators:
#   - TWFE (naive pooled)
#   - Gardner (2022) two-stage
#   - L0-penalised greedy (10 penalty values)
#   - DP Bayes (collapsed Gibbs, diffuse priors)
#
# Output:
#   - Bias and variance for each estimator on the overall ATT and CATTs
#   - Single summary table written to results/sim_summary.csv
# ============================================================

library(fixest)
library(data.table)

# ---------- Constants ----------

N_MC          <- 100
T_PERIODS     <- 10
TREATED_GS    <- c(3, 5, 7)
SIGMA_EPS     <- 1.0
STRUCTURES    <- c("homog", "six_partition", "full_het")

# Cohort-size specifications. Each entry is a named vector of length 4 giving
# the number of units assigned to (never-treated, cohort 3, cohort 5, cohort 7).
# uniform_*  -> all four cohorts the same size;
# unequal_*  -> heterogeneous sizes across cohorts.
COHORT_SIZE_SPECS <- list(
  uniform_small = c("0" = 10, "3" = 10, "5" = 10, "7" = 10),
  unequal_small = c("0" =  5, "3" = 10, "5" = 13, "7" = 15)
)

# DP Bayes
N_GIBBS       <- 500
N_BURNIN      <- 100              # first 20% of iterations discarded as burn-in
ALPHA_DP_GRID <- c(1.0, 20.0)     # concentration: 1 favours few clusters, 20 favours many
SIGMA0_SQ     <- 100.0            # diffuse base measure
MU0           <- 0.0

# L0 penalty grid: 10 log-spaced values spanning extreme cases
L_GRID <- 10^seq(-2, 2, length.out = 10)

set.seed(42)
MC_SEEDS <- sample.int(1e7, N_MC)

# ---------- CATT cells ----------

build_catt_cells <- function() {
  cells <- expand.grid(cohort = TREATED_GS, time = seq_len(T_PERIODS))
  cells <- cells[cells$time >= cells$cohort, ]
  cells <- cells[order(cells$cohort, cells$time), ]
  cells$cell_id <- seq_len(nrow(cells))
  rownames(cells) <- NULL
  cells
}

CATT_CELLS <- build_catt_cells()
K          <- nrow(CATT_CELLS)   # 8 + 6 + 4 = 18

# ---------- Calendar-time ATT aggregation (Callaway & Sant'Anna 2021) ----------
#
# Replaces the simple equal-weight mean with the heterogeneity-preserving
# calendar-time aggregator theta^O_c, eqs. (3.8) and (3.12) of the paper:
#
#   theta_c(t)  = sum_g 1{t >= g} P(G=g | G <= t) * ATT(g, t)
#   theta^O_c   = (1 / (T-1)) * sum_{t=2}^{T} theta_c(t)
#
# Cohort sizes are uniform within each scenario in this simulation, so
# P(G=g | G <= t) reduces to 1 / #{treated cohorts with g <= t}.

aggregate_calendar_time <- function(catt_vec, cohort_sizes = NULL) {
  if (is.null(cohort_sizes)) cohort_sizes <- rep(1, length(TREATED_GS))
  treated_n <- setNames(cohort_sizes, as.character(TREATED_GS))
  theta_c <- numeric(T_PERIODS)
  for (t in seq_len(T_PERIODS)) {
    g_treated <- TREATED_GS[TREATED_GS <= t]
    if (length(g_treated) == 0L) { theta_c[t] <- 0; next }
    w <- treated_n[as.character(g_treated)]
    w <- w / sum(w)
    att_gt <- vapply(g_treated, function(g) {
      idx <- which(CATT_CELLS$cohort == g & CATT_CELLS$time == t)
      if (length(idx) == 0L) NA_real_ else catt_vec[idx]
    }, numeric(1))
    theta_c[t] <- sum(w * att_gt, na.rm = TRUE)
  }
  sum(theta_c[2:T_PERIODS], na.rm = TRUE) / (T_PERIODS - 1)
}

# Event-time aggregation, CS eqs. (3.4) and (3.12):
#   theta_es(e) = sum_g 1{g+e <= T} P(G=g | G+e <= T) ATT(g, g+e)
#   theta^O_es  = (1 / (T-1)) sum_{e=0}^{T-2} theta_es(e)
# Uniform cohort sizes => P(G=g | G+e <= T) = 1 / #{eligible cohorts at e}.

aggregate_event_time <- function(catt_vec, cohort_sizes = NULL) {
  if (is.null(cohort_sizes)) cohort_sizes <- rep(1, length(TREATED_GS))
  treated_n <- setNames(cohort_sizes, as.character(TREATED_GS))
  theta_es <- numeric(T_PERIODS - 1)        # event times e = 0,...,T-2
  for (e_idx in seq_along(theta_es)) {
    e <- e_idx - 1L
    g_eligible <- TREATED_GS[TREATED_GS + e <= T_PERIODS]
    if (length(g_eligible) == 0L) { theta_es[e_idx] <- 0; next }
    w <- treated_n[as.character(g_eligible)]
    w <- w / sum(w)
    att_ge <- vapply(g_eligible, function(g) {
      idx <- which(CATT_CELLS$cohort == g & CATT_CELLS$time == g + e)
      if (length(idx) == 0L) NA_real_ else catt_vec[idx]
    }, numeric(1))
    theta_es[e_idx] <- sum(w * att_ge, na.rm = TRUE)
  }
  sum(theta_es, na.rm = TRUE) / (T_PERIODS - 1)
}

# Selection-style overall ATT, CS eqs. (3.7) and (3.11):
#   theta_sel(g) = (1 / (T - g + 1)) sum_{t=g}^T ATT(g, t)
#   theta^O_sel  = sum_g theta_sel(g) P(G=g | G <= T)
# Uniform cohort sizes => weights are uniform across treated cohorts.

aggregate_selection <- function(catt_vec, cohort_sizes = NULL) {
  if (is.null(cohort_sizes)) cohort_sizes <- rep(1, length(TREATED_GS))
  theta_sel <- vapply(TREATED_GS, function(g) {
    cells_g <- which(CATT_CELLS$cohort == g)
    mean(catt_vec[cells_g], na.rm = TRUE)
  }, numeric(1))
  w <- cohort_sizes / sum(cohort_sizes)
  sum(w * theta_sel, na.rm = TRUE)
}

# ---------- True CATT vector by scenario ----------

make_true_catts <- function(structure) {
  if (structure == "homog") {
    catts    <- rep(2, K)
    group_id <- rep(1L, K)
  } else if (structure == "six_partition") {
    # Mixed-gap design:
    #   trio (close):       1.0, 1.2, 1.4
    #   duo  (medium-far):  4.0, 4.2
    #   singleton (far):    8.0
    group_vals <- c(1.0, 1.2, 1.4, 4.0, 4.2, 8.0)
    group_id   <- rep(1:6, each = 3)        # K = 18 = 6 * 3
    catts      <- group_vals[group_id]
  } else if (structure == "full_het") {
    catts    <- seq(1, 8, length.out = K)
    group_id <- seq_len(K)
  } else stop("Unknown structure: ", structure)

  list(catts = catts, group_id = group_id)
}

# ---------- DGP ----------

generate_panel <- function(cohort_sizes_vec, true_catts, seed) {
  # cohort_sizes_vec: named numeric of length 4 with names matching the cohort
  # indicators used elsewhere, i.e. "0", "3", "5", "7".
  set.seed(seed)

  cohorts_all <- c(0, TREATED_GS)
  sizes_in_order <- cohort_sizes_vec[as.character(cohorts_all)]
  N_total <- sum(sizes_in_order)

  unit_ids    <- seq_len(N_total)
  unit_cohort <- rep(cohorts_all, times = sizes_in_order)
  alpha_i     <- rnorm(N_total, mean = 0, sd = 1)

  panel <- CJ(unit = unit_ids, time = seq_len(T_PERIODS))
  panel[, cohort := unit_cohort[unit]]
  panel[, alpha  := alpha_i[unit]]
  panel[, lambda := 0.1 * (time - 1)]
  panel[, D      := as.integer(cohort > 0 & time >= cohort)]

  catt_lookup <- as.data.table(CATT_CELLS)
  catt_lookup[, tau_true := true_catts]

  panel <- merge(panel, catt_lookup[, .(cohort, time, tau_true)],
                 by = c("cohort", "time"), all.x = TRUE)
  panel[is.na(tau_true), tau_true := 0]

  panel[, eps := rnorm(.N, mean = 0, sd = SIGMA_EPS)]
  panel[, y   := alpha + lambda + D * tau_true + eps]

  setorder(panel, unit, time)

  for (k in seq_len(K)) {
    g_k <- CATT_CELLS$cohort[k]
    t_k <- CATT_CELLS$time[k]
    panel[, paste0("D_", k) := as.integer(cohort == g_k & time == t_k)]
  }

  panel[]
}

# ---------- Within transformation ----------

within_transform <- function(panel) {
  cols <- c("y", paste0("D_", seq_len(K)))
  out  <- copy(panel)
  for (col in cols) {
    out[, paste0(col, "_unit")  := mean(get(col)), by = unit]
    out[, paste0(col, "_time")  := mean(get(col)), by = time]
    grand <- mean(out[[col]])
    out[, paste0(col, "_dm")    := get(col) - get(paste0(col, "_unit")) -
                                    get(paste0(col, "_time")) + grand]
  }
  out
}

# ---------- Estimator: TWFE (naive pooled) ----------

estimate_twfe <- function(panel) {
  fit     <- feols(y ~ D | unit + time, data = panel)
  tau_hat <- as.numeric(coef(fit)["D"])
  catt    <- rep(tau_hat, K)
  list(catt = catt, n_partitions = 1L)
}

# ---------- Estimator: Fully flexible (used by L0 and DP) ----------

estimate_flexible <- function(panel_dm) {
  X <- as.matrix(panel_dm[, paste0("D_", seq_len(K), "_dm"), with = FALSE])
  y <- panel_dm$y_dm
  XtX_inv <- solve(crossprod(X))
  beta    <- as.numeric(XtX_inv %*% crossprod(X, y))
  resid   <- y - X %*% beta
  rss     <- sum(resid^2)
  list(catt = beta, rss = rss, XtX_inv = XtX_inv, X = X, y = y,
       n_partitions = K)
}

# ---------- Estimator: Gardner (2022) two-stage ----------

estimate_gardner <- function(panel) {
  fit_y0 <- feols(y ~ 1 | unit + time, data = panel[D == 0])
  panel_local <- copy(panel)
  panel_local[, y0_hat := predict(fit_y0, newdata = panel_local)]
  panel_local[, resid  := y - y0_hat]

  catt_df <- panel_local[D == 1, .(catt_est = mean(resid)), by = .(cohort, time)]
  cells   <- merge(as.data.table(CATT_CELLS), catt_df,
                   by = c("cohort", "time"), all.x = TRUE)
  setorder(cells, cell_id)

  list(catt = cells$catt_est,
       n_partitions = sum(!is.na(cells$catt_est)))
}

# ---------- Estimator: L0 greedy (multiple L values) ----------

estimate_l0_grid <- function(flex_catts, n_per_cell, L_values) {
  results <- vector("list", length(L_values))
  names(results) <- paste0("L_", formatC(L_values, format = "g", digits = 4))

  for (idx in seq_along(L_values)) {
    L      <- L_values[idx]
    groups <- as.list(seq_len(K))
    n_grp  <- n_per_cell
    tau_grp <- flex_catts

    repeat {
      m         <- length(groups)
      if (m == 1L) break
      best_dobj <- Inf
      best_pair <- NULL

      for (i in seq_len(m - 1)) {
        for (j in (i + 1):m) {
          cross_pairs <- length(groups[[i]]) * length(groups[[j]])
          ni <- n_grp[i]; nj <- n_grp[j]
          delta_rss   <- (ni * nj / (ni + nj)) * (tau_grp[i] - tau_grp[j])^2
          delta_obj   <- delta_rss - L * cross_pairs
          if (delta_obj < best_dobj) {
            best_dobj <- delta_obj
            best_pair <- c(i, j)
          }
        }
      }
      if (best_dobj >= 0) break

      i <- best_pair[1]; j <- best_pair[2]
      ni <- n_grp[i]; nj <- n_grp[j]
      pooled <- (ni * tau_grp[i] + nj * tau_grp[j]) / (ni + nj)

      groups[[i]]  <- c(groups[[i]], groups[[j]])
      tau_grp[i]   <- pooled
      n_grp[i]     <- ni + nj
      groups[[j]]  <- NULL
      tau_grp      <- tau_grp[-j]
      n_grp        <- n_grp[-j]
    }

    catt_est <- numeric(K)
    for (g in seq_along(groups)) catt_est[groups[[g]]] <- tau_grp[g]

    results[[idx]] <- list(
      L        = L,
      catt     = catt_est,
      n_groups = length(groups)
    )
  }
  results
}

# ---------- Estimator: DP Bayes (collapsed Gibbs) ----------

estimate_dp_bayes <- function(panel_dm, sigma2_hat, alpha) {
  X  <- as.matrix(panel_dm[, paste0("D_", seq_len(K), "_dm"), with = FALSE])
  y  <- panel_dm$y_dm
  NT <- length(y)
  kappa <- SIGMA0_SQ / sigma2_hat

  log_marg_lik <- function(z_relab) {
    m <- max(z_relab)
    Xg <- matrix(0, nrow = NT, ncol = m)
    for (l in seq_len(m)) Xg[, l] <- rowSums(X[, z_relab == l, drop = FALSE])
    A   <- diag(m) + kappa * crossprod(Xg)
    Xty <- crossprod(Xg, y)
    -0.5 * as.numeric(determinant(A, logarithm = TRUE)$modulus) +
      kappa / (2 * sigma2_hat) * as.numeric(t(Xty) %*% solve(A, Xty))
  }

  z <- seq_len(K)                                        # singletons
  theta_draws       <- matrix(NA_real_, nrow = N_GIBBS, ncol = K)
  n_clusters_trace  <- integer(N_GIBBS)

  for (iter in seq_len(N_GIBBS)) {

    for (k in seq_len(K)) {
      z_minus_k       <- z[-k]
      unique_clusters <- sort(unique(z_minus_k))
      n_minus_k       <- tabulate(z_minus_k, nbins = max(z))
      c_new           <- max(z) + 1L
      candidates      <- c(unique_clusters, c_new)
      log_probs       <- numeric(length(candidates))

      for (idx in seq_along(candidates)) {
        c_cand  <- candidates[idx]
        z_trial <- z; z_trial[k] <- c_cand
        z_trial <- match(z_trial, sort(unique(z_trial)))
        crp_wt  <- if (c_cand == c_new) alpha else n_minus_k[c_cand]
        log_probs[idx] <- log(crp_wt) + log_marg_lik(z_trial)
      }

      log_probs <- log_probs - max(log_probs)
      probs     <- exp(log_probs) / sum(exp(log_probs))
      z[k]      <- candidates[sample(length(candidates), 1L, prob = probs)]
      z         <- match(z, sort(unique(z)))
    }

    m  <- max(z)
    n_clusters_trace[iter] <- m
    Xg <- matrix(0, nrow = NT, ncol = m)
    for (l in seq_len(m)) Xg[, l] <- rowSums(X[, z == l, drop = FALSE])
    Sigma_star <- solve(crossprod(Xg) / sigma2_hat + diag(m) / SIGMA0_SQ)
    mu_star    <- Sigma_star %*% (crossprod(Xg, y) / sigma2_hat +
                                  MU0 * rep(1, m) / SIGMA0_SQ)
    L_chol <- tryCatch(chol(Sigma_star),
                       error = function(e) chol(Sigma_star + 1e-8 * diag(m)))
    phi <- as.numeric(mu_star) + as.numeric(t(L_chol) %*% rnorm(m))
    theta_draws[iter, ] <- phi[z]
  }

  keep <- (N_BURNIN + 1L):N_GIBBS
  catt_mean <- colMeans(theta_draws[keep, , drop = FALSE])
  list(catt = catt_mean,
       n_partitions = mean(n_clusters_trace[keep]))
}

# ---------- Main MC loop ----------

scenario_grid <- expand.grid(
  cohort_spec = names(COHORT_SIZE_SPECS),
  structure   = STRUCTURES,
  stringsAsFactors = FALSE
)
# Homogeneous truth is invariant to cohort sizes, so we only run it for
# the uniform spec and drop the redundant unequal-cohort homogeneous row.
scenario_grid <- subset(scenario_grid,
                        !(structure == "homog" & cohort_spec != "uniform_small"))
rownames(scenario_grid) <- NULL

results_all <- list()
row_counter <- 1L

for (s in seq_len(nrow(scenario_grid))) {

  cohort_spec <- scenario_grid$cohort_spec[s]
  structure   <- scenario_grid$structure[s]
  sizes_vec   <- COHORT_SIZE_SPECS[[cohort_spec]]              # length 4: never + 3 treated
  treated_n   <- sizes_vec[as.character(TREATED_GS)]            # length 3: treated only

  truth        <- make_true_catts(structure)
  true_catts   <- truth$catts
  true_att_c   <- aggregate_calendar_time(true_catts, treated_n)
  true_att_es  <- aggregate_event_time(true_catts,    treated_n)
  true_att_sel <- aggregate_selection(true_catts,     treated_n)
  true_n_partitions <- length(unique(truth$group_id))
  n_per_cell <- vapply(seq_len(K), function(k)
    treated_n[as.character(CATT_CELLS$cohort[k])], numeric(1))

  cat(sprintf("Scenario %d/%d: cohort_spec=%s, structure=%s\n",
              s, nrow(scenario_grid), cohort_spec, structure))

  rep_storage <- list()

  for (rep in seq_len(N_MC)) {

    panel    <- generate_panel(sizes_vec, true_catts, MC_SEEDS[rep])
    panel_dm <- within_transform(panel)

    flex     <- estimate_flexible(panel_dm)
    df_resid <- nrow(panel_dm) - length(unique(panel$unit)) -
                length(unique(panel$time)) + 1 - K
    sigma2_hat <- flex$rss / df_resid

    twfe_out    <- estimate_twfe(panel)
    gardner_out <- estimate_gardner(panel)
    l0_out      <- estimate_l0_grid(flex$catt, n_per_cell, L_GRID)
    dp_out_list <- lapply(ALPHA_DP_GRID, function(a)
      estimate_dp_bayes(panel_dm, sigma2_hat, alpha = a))
    names(dp_out_list) <- sprintf("alpha=%g", ALPHA_DP_GRID)

    rep_row <- data.table()
    rep_row <- rbind(rep_row, data.table(
      method = "TWFE", catt = list(twfe_out$catt),
      n_partitions = twfe_out$n_partitions))
    rep_row <- rbind(rep_row, data.table(
      method = "Flexible", catt = list(flex$catt),
      n_partitions = flex$n_partitions))
    rep_row <- rbind(rep_row, data.table(
      method = "Gardner", catt = list(gardner_out$catt),
      n_partitions = gardner_out$n_partitions))
    for (a_idx in seq_along(ALPHA_DP_GRID)) {
      a_val <- ALPHA_DP_GRID[a_idx]
      rep_row <- rbind(rep_row, data.table(
        method = sprintf("DP Bayes (alpha=%g)", a_val),
        catt   = list(dp_out_list[[a_idx]]$catt),
        n_partitions = dp_out_list[[a_idx]]$n_partitions))
    }
    for (idx in seq_along(L_GRID)) {
      rep_row <- rbind(rep_row, data.table(
        method = sprintf("L0 (L=%.4g)", L_GRID[idx]),
        catt   = list(l0_out[[idx]]$catt),
        n_partitions = l0_out[[idx]]$n_groups))
    }

    rep_row[, replication := rep]
    rep_row[, cohort_spec := cohort_spec]
    rep_row[, structure   := structure]

    rep_storage[[rep]] <- rep_row
  }

  results_all[[row_counter]] <- list(
    scenario          = scenario_grid[s, ],
    cohort_spec       = cohort_spec,
    treated_n         = treated_n,
    true_catts        = true_catts,
    true_att_c        = true_att_c,
    true_att_es       = true_att_es,
    true_att_sel      = true_att_sel,
    true_n_partitions = true_n_partitions,
    raw               = rbindlist(rep_storage)
  )
  row_counter <- row_counter + 1L
}

# ---------- Aggregate: bias and variance ----------

summary_rows <- list()

for (item in results_all) {

  scen         <- item$scenario
  treated_n    <- item$treated_n
  true_catts   <- item$true_catts
  true_att_c   <- item$true_att_c
  true_att_es  <- item$true_att_es
  true_att_sel <- item$true_att_sel
  true_m       <- item$true_n_partitions
  raw          <- item$raw

  for (mth in unique(raw$method)) {
    sub <- raw[method == mth]

    # Recompute the three ATTs from each replication's CATT vector,
    # using the cohort sizes that match the scenario.
    att_c_vec   <- vapply(sub$catt, aggregate_calendar_time, numeric(1),
                          cohort_sizes = treated_n)
    att_es_vec  <- vapply(sub$catt, aggregate_event_time,    numeric(1),
                          cohort_sizes = treated_n)
    att_sel_vec <- vapply(sub$catt, aggregate_selection,     numeric(1),
                          cohort_sizes = treated_n)

    att_c_bias    <- mean(att_c_vec,   na.rm = TRUE) - true_att_c
    att_c_var     <- var(att_c_vec,    na.rm = TRUE)
    att_es_bias   <- mean(att_es_vec,  na.rm = TRUE) - true_att_es
    att_es_var    <- var(att_es_vec,   na.rm = TRUE)
    att_sel_bias  <- mean(att_sel_vec, na.rm = TRUE) - true_att_sel
    att_sel_var   <- var(att_sel_vec,  na.rm = TRUE)

    catt_mat  <- do.call(rbind, sub$catt)        # N_MC x K
    catt_mean <- colMeans(catt_mat, na.rm = TRUE)
    catt_var  <- apply(catt_mat, 2, var, na.rm = TRUE)
    catt_bias <- catt_mean - true_catts

    avg_abs_bias_catt <- mean(abs(catt_bias), na.rm = TRUE)
    avg_var_catt      <- mean(catt_var, na.rm = TRUE)
    rmse_catt         <- sqrt(mean(catt_bias^2 + catt_var, na.rm = TRUE))
    avg_n_partitions  <- mean(sub$n_partitions, na.rm = TRUE)

    summary_rows[[length(summary_rows) + 1]] <- data.table(
      cohort_spec        = scen$cohort_spec,
      structure          = scen$structure,
      true_partitions    = true_m,
      method             = mth,
      avg_n_partitions   = avg_n_partitions,
      ATT_c_bias         = att_c_bias,
      ATT_c_variance     = att_c_var,
      ATT_es_bias        = att_es_bias,
      ATT_es_variance    = att_es_var,
      ATT_sel_bias       = att_sel_bias,
      ATT_sel_variance   = att_sel_var,
      CATT_avg_abs_bias  = avg_abs_bias_catt,
      CATT_avg_variance  = avg_var_catt,
      CATT_rmse          = rmse_catt
    )
  }
}

summary_df <- rbindlist(summary_rows)

method_levels <- c("TWFE", "Gardner", "Flexible",
                   sprintf("DP Bayes (alpha=%g)", ALPHA_DP_GRID),
                   sprintf("L0 (L=%.4g)", L_GRID))
summary_df[, method := factor(method, levels = method_levels)]
summary_df[, cohort_spec := factor(cohort_spec, levels = names(COHORT_SIZE_SPECS))]
setorder(summary_df, cohort_spec, structure, method)

num_cols <- c("avg_n_partitions",
              "ATT_c_bias", "ATT_c_variance",
              "ATT_es_bias", "ATT_es_variance",
              "ATT_sel_bias", "ATT_sel_variance",
              "CATT_avg_abs_bias", "CATT_avg_variance", "CATT_rmse")
summary_df[, (num_cols) := lapply(.SD, round, 4), .SDcols = num_cols]

# ---------- Final summary table ----------

print(summary_df)

dir.create("results", showWarnings = FALSE)
fwrite(summary_df, "results/sim_summary.csv")

cat("\nSaved summary to results/sim_summary.csv\n")
cat("Scenarios run:", nrow(scenario_grid), "\n")
cat("Methods compared:", length(method_levels),
    "(including", length(L_GRID), "L0 penalty values)\n")
cat("Replications per scenario:", N_MC, "\n")
