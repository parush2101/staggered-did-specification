# staggered-did-specification

**ℓ₀-Penalized CATT Estimation for Staggered Difference-in-Differences**

---

## Motivation

In a staggered DiD design, units enter treatment at different calendar times, generating a set of **Cohort-Average Treatment effects on the Treated (CATTs)** — one for each (cohort, time) pair. The fully flexible estimator treats every CATT as a distinct parameter, which is unbiased but inefficient when some CATTs are equal. The fully pooled (naive TWFE) estimator imposes a single treatment effect across all cohorts and periods, which is efficient but biased when treatment effects are genuinely heterogeneous.

This project develops a **data-driven model specification approach**: rather than committing to either extreme, we penalise the *number of distinct CATT values* in the model, allowing the data to determine which CATTs can be pooled.

---

## Method

### Setup

Consider a balanced panel with $N$ units and $T$ time periods. Units belong to one of $G$ treated cohorts (entering treatment at period $g$) or a never-treated group. Define the cohort-time treatment dummy:

$$D_{git} = \mathbf{1}\{\text{unit } i \in \text{cohort } g,\; t = s\}$$

The fully flexible model estimates one CATT per (cohort, calendar-time) cell:

$$y_{it} = \sum_{g} \sum_{t \geq g} \tau_{gt} \cdot D_{git} + \alpha_i + \lambda_t + \varepsilon_{it}$$

### ℓ₀ Penalized Objective

A **partition** $P = \{G_1, \ldots, G_m\}$ of all $K$ CATTs assigns each CATT to a group whose members share a common parameter. The penalized objective is:

$$Q(P) = \text{RSS}(P) + L \cdot \underbrace{{(j,k) : j \text{ and } k \text{ in different groups}\}}_{\text{number of cross-group CATT pairs}}$$

- When $L = 0$: fully flexible (minimise RSS, no penalty for having $K$ distinct parameters).
- When $L \to \infty$: fully pooled (penalty dominates, all CATTs collapsed to one group).
- For intermediate $L$: groups of equal CATTs are recovered from the data.

### Greedy Algorithm

Exact minimisation of $Q(P)$ is NP-hard (the search space is the set of all partitions of $K$ elements, the Bell number $B(K)$). We use a **greedy merging** heuristic:

1. Initialise with $K$ singleton groups (fully flexible).
2. For every pair of current groups $(A, B)$, compute the net change in objective if they are merged:
   $$\Delta \text{Obj}(A,B) = \underbrace{\Delta \text{RSS}_{A \cup B}}_{\geq\, 0} - \underbrace{L \cdot |A| \cdot |B|}_{\text{penalty saving}}$$
3. Accept the merge with the most-negative $\Delta \text{Obj}$ (if $< 0$).
4. Repeat until no beneficial merge remains.

**Complexity:** $O(K^2)$ objective evaluations per step, $O(K)$ steps → $O(K^3)$ total. Each objective evaluation is a small OLS after a one-time within-transformation (TWFE double-demeaning), making this fast in practice.

**Penalty interpretation:** Two CATTs $\tau_j$ and $\tau_k$ are merged only if the RSS increase from constraining them equal is smaller than $L$. For singleton groups, this threshold is exactly $L$.

---

## Simulation Design

| Parameter | Value |
|---|---|
| Units per cohort | 100 |
| Time periods | 8 |
| Cohorts | 3, 5, 7 (+ never-treated) |
| Unit FE | $\alpha_i \sim \mathcal{N}(0,1)$ |
| Time FE | $\lambda_t = 0.1(t-1)$ |
| Noise | $\varepsilon_{it} \sim \mathcal{N}(0, 0.25)$ |

**True CATT structure (2 groups):**

| Group | CATTs | True $\tau$ |
|---|---|---|
| A | Cohort 3 at $t=3,\ldots,8$ and Cohort 5 at $t=5,\ldots,8$ | 2 |
| B | Cohort 7 at $t=7, 8$ | 4 |

Total: 12 CATTs, true partition has 2 groups.

**Estimators compared:**
- **Flexible** — 12 parameters, no pooling
- **L0 Greedy** — penalized, for $L \in \{0.5,\, 2,\, 10,\, 50\}$
- **Pooled** — 1 parameter, full pooling

---

## Repository Structure

```
staggered-did-specification/
├── README.md
├── .gitignore
├── R/
│   ├── dgp.R           # Data generating process and treatment dummies
│   ├── estimators.R    # Flexible, pooled, and greedy L0 estimators
│   └── plots.R         # ggplot2 visualisation functions
└── simulations/
    └── run_simulation.R  # Main script: generates data, runs all estimators,
                          #              prints summary table, saves plots
```

---

## How to Run

```r
# From the project root directory
source("simulations/run_simulation.R")
```

Output:
- Console summary table (# groups, RSS, RMSE per method)
- `simulations/catt_estimates.png` — CATT estimates by cohort
- `simulations/solution_path.png` — # groups vs. penalty L

**Requirements:** R ≥ 4.0, packages `ggplot2`, `scales`, `MASS` (base).

---

## Key Result

For this DGP, the greedy algorithm recovers the true 2-group partition for a wide range of $L$ values because:

- The RSS increase from *correctly* merging two equal-effect CATTs is small ($\approx \sigma^2 / N$).
- The RSS increase from *incorrectly* merging CATTs with different true effects ($\tau=2$ vs. $\tau=4$) is large ($\approx N \cdot (\Delta\tau)^2$).

As long as $L$ falls between these two thresholds, the algorithm terminates at the true partition.

---

## References

- Callaway, B. & Sant'Anna, P. H. C. (2021). Difference-in-differences with multiple time periods. *Journal of Econometrics*, 225, 200–230.
- Sun, L. & Abraham, S. (2021). Estimating dynamic treatment effects in event studies with heterogeneous treatment effects. *Journal of Econometrics*, 225, 175–199.
- Wooldridge, J. M. (2025). Two-way fixed effects, the two-way Mundlak regression, and difference-in-differences estimators. *Empirical Economics*, 69, 2545–2587.
