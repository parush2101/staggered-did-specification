# =============================================================================
# estimators.R — CATT Estimators
# =============================================================================
#
# Three estimators are implemented:
#
#   1. estimate_flexible()  — fully flexible: each (g,t) cell is a separate
#                             parameter. Equivalent to Sun & Abraham (2021)
#                             / Callaway & Sant'Anna (2021).
#
#   2. estimate_pooled()    — fully pooled: one common ATT across all cohorts
#                             and time periods. Equivalent to naive TWFE.
#
#   3. greedy_l0()          — ℓ₀-penalized greedy estimator. Starts from the
#                             fully flexible partition and iteratively merges
#                             pairs of CATT groups, accepting a merge only when
#                             it reduces the penalized objective:
#
#                             Q(P) = RSS(P) + L · #{pairs (j,k) : j,k in diff groups}
#
#                             A merge of groups A (size nA) and B (size nB) is
#                             accepted iff  ΔRSS < L · nA · nB.
#
# All estimators use the TWFE within-transformation (double-demeaning) for
# speed.  The Frisch–Waugh–Lovell theorem guarantees this gives identical
# point estimates to including unit and time dummies explicitly.
# =============================================================================

source("R/dgp.R")

# =============================================================================
# Within transformation helpers
# =============================================================================

# Double-demean a single variable for a balanced panel (one-step exact).
#   x_dm = x - unit_mean(x) - time_mean(x) + grand_mean(x)
double_demean <- function(x, unit, time) {
  gm <- mean(x)
  x - ave(x, unit) - ave(x, time) + gm
}


# Pre-compute the within-transformed outcome and all CATT dummies once.
# Stores results as y_dm and D_<g>_<t>_dm columns.
prepare_within <- function(panel) {
  panel$y_dm <- double_demean(panel$y, panel$unit, panel$time)

  catt_cols <- grep("^D_", names(panel), value = TRUE)
  for (col in catt_cols) {
    panel[[paste0(col, "_dm")]] <- double_demean(panel[[col]], panel$unit, panel$time)
  }
  return(panel)
}


# =============================================================================
# Core estimator: OLS for a given partition
# =============================================================================

# Given a partition (list of character vectors of CATT labels), build the
# group-level regressors by summing the within-transformed CATT dummies within
# each group, then run OLS of y_dm on those group regressors.
#
# Returns: list(coefs, rss, partition)
#   coefs     — named numeric vector, one coefficient per group
#   rss       — residual sum of squares
#   partition — the partition that was estimated
estimate_partition <- function(panel_dm, partition) {

  # Build group regressors (sum of within-transformed dummies in each group)
  X_list <- lapply(seq_along(partition), function(k) {
    dm_cols <- paste0(partition[[k]], "_dm")
    rowSums(panel_dm[, dm_cols, drop = FALSE])
  })
  X <- do.call(cbind, X_list)                 # N*T × #groups
  y <- panel_dm$y_dm                          # N*T × 1

  # OLS via normal equations (small system: #groups << N*T)
  XtX  <- crossprod(X)
  Xty  <- crossprod(X, y)
  coefs <- tryCatch(
    solve(XtX, Xty),
    error = function(e) MASS::ginv(XtX) %*% Xty  # fallback for near-singular
  )

  resids <- y - X %*% coefs
  rss    <- sum(resids^2)

  names(coefs) <- paste0("G_", seq_along(partition))
  return(list(coefs = coefs, rss = rss, partition = partition))
}


# =============================================================================
# Named estimators
# =============================================================================

# Fully flexible: each CATT is its own group
estimate_flexible <- function(panel_dm) {
  catts     <- get_catt_list(max(panel_dm$time))
  partition <- as.list(catts$label)           # K singleton groups
  estimate_partition(panel_dm, partition)
}

# Fully pooled: all CATTs in one group
estimate_pooled <- function(panel_dm) {
  catts     <- get_catt_list(max(panel_dm$time))
  partition <- list(catts$label)              # 1 group containing all CATTs
  estimate_partition(panel_dm, partition)
}


# =============================================================================
# ℓ₀ greedy estimator
# =============================================================================

# greedy_l0()
#
# Algorithm:
#   1. Initialise: one singleton group per CATT (fully flexible).
#   2. For every pair of current groups (A, B):
#        Compute ΔRSS   = RSS(merged partition) − RSS(current partition)
#        Compute saving = L · |A| · |B|          (penalty reduction)
#        Net change     = ΔRSS − saving
#   3. Accept the merge with the most-negative net change (if < 0).
#   4. Repeat until no beneficial merge exists.
#
# Complexity: O(K² · lm_cost) per iteration, O(K) iterations → O(K³) total.
# For K = 12 this is very fast.
#
# Arguments:
#   panel_dm  — panel with within-transformed columns (from prepare_within)
#   L         — penalty per CATT pair in different groups
#   verbose   — print merge log
greedy_l0 <- function(panel_dm, L, verbose = TRUE) {

  catts       <- get_catt_list(max(panel_dm$time))
  K           <- nrow(catts)
  partition   <- as.list(catts$label)         # start: K singletons
  group_sizes <- rep(1L, K)                   # number of CATTs per group

  current_rss <- estimate_partition(panel_dm, partition)$rss

  if (verbose) {
    cat(sprintf("L = %.2f | Start: %d groups | RSS = %.4f\n",
                L, K, current_rss))
  }

  step    <- 0L
  changed <- TRUE

  while (changed && length(partition) > 1L) {
    changed   <- FALSE
    best_dobj <- 0          # only accept strictly negative changes
    best_move <- NULL
    n_grp     <- length(partition)

    for (a in seq_len(n_grp - 1L)) {
      for (b in seq(a + 1L, n_grp)) {

        # Proposed partition with groups a and b merged
        trial            <- partition
        trial[[a]]       <- c(trial[[a]], trial[[b]])
        trial[[b]]       <- NULL

        trial_rss        <- estimate_partition(panel_dm, trial)$rss
        delta_rss        <- trial_rss - current_rss      # >= 0

        penalty_saving   <- L * group_sizes[a] * group_sizes[b]
        delta_obj        <- delta_rss - penalty_saving   # <0 → merge beneficial

        if (delta_obj < best_dobj) {
          best_dobj <- delta_obj
          best_move <- list(
            a         = a,
            b         = b,
            partition = trial,
            rss       = trial_rss,
            sizes     = c(group_sizes[-c(a, b)],
                          group_sizes[a] + group_sizes[b])
          )
        }
      }
    }

    if (!is.null(best_move)) {
      step        <- step + 1L
      partition   <- best_move$partition
      current_rss <- best_move$rss
      group_sizes <- best_move$sizes
      changed     <- TRUE

      if (verbose) {
        cat(sprintf("  Step %d: merged groups %d & %d → %d groups | RSS = %.4f | ΔObj = %.4f\n",
                    step, best_move$a, best_move$b,
                    length(partition), current_rss, best_dobj))
      }
    }
  }

  final <- estimate_partition(panel_dm, partition)

  if (verbose) {
    cat(sprintf("  Done: %d groups | RSS = %.4f\n\n",
                length(partition), current_rss))
  }

  return(list(
    partition = partition,
    n_groups  = length(partition),
    rss       = current_rss,
    coefs     = final$coefs,
    L         = L,
    n_steps   = step
  ))
}


# =============================================================================
# Helper: expand partition results to CATT-level data frame
# =============================================================================

# Returns a data frame with one row per CATT, with columns:
#   cohort, time, label, tau_true, estimate, group, method
collect_estimates <- function(result, panel_dm, method_label) {
  catts     <- get_catt_list(max(panel_dm$time))
  partition <- result$partition
  coefs     <- result$coefs

  out <- catts
  out$estimate <- NA_real_
  out$group    <- NA_integer_
  out$method   <- method_label

  for (k in seq_along(partition)) {
    idx             <- which(out$label %in% partition[[k]])
    out$estimate[idx] <- unname(coefs[k])
    out$group[idx]  <- k
  }
  return(out)
}
