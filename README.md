# staggered-did-specification

**Specification and Partition Selection for Staggered Difference-in-Differences**

A Bayesian (Dirichlet-Process) approach to deciding which cohort-time treatment
effects (CATTs) to pool in a staggered DiD design. The paper is
`staggered_did_specification.tex` (compiled: `staggered_did_specification.pdf`).

---

## Motivation

In a staggered DiD design, units enter treatment at different calendar times,
generating a set of Cohort-Average Treatment effects on the Treated (CATTs), one
per (cohort, time) cell. Estimating every CATT separately (the fully flexible
estimator) is unbiased but inefficient when some CATTs are equal; pooling all of
them into a single TWFE coefficient is efficient but biased when the effects are
genuinely heterogeneous. This project frames the choice between these extremes as
a partition-selection problem and solves it with a Dirichlet-Process (DP) mixture
prior over the CATTs.

---

## Method

The primary estimator is a **full Bayesian Dirichlet-Process mixture** over the
CATTs. A collapsed Gibbs sampler (drawing the error variance from its
inverse-gamma full conditional each sweep) delivers:

- posterior point estimates,
- credible intervals that marginalize the unknown partition, and
- co-clustering probabilities for every pair of CATTs.

Honest inference under partition selection comes from averaging over partitions
rather than conditioning on a single grouping.

**ℓ₀ connection (remark).** Holding the error variance fixed, the model's maximum
a posteriori partition coincides with an ℓ₀-penalized least-squares estimator (a
BIC-type criterion when the variance is free), linking the Bayesian formulation
to frequentist homogeneity pursuit. This estimator is carried through the
simulations as a comparison and behaves much like the posterior.

---

## Repository structure

```
staggered-did-specification/
├── staggered_did_specification.tex   # the paper
├── staggered_did_specification.pdf   # compiled paper
├── references.bib
├── staggered_sim.R                   # Monte Carlo study (Parts 1-3)
├── mpdta_run.R                       # Callaway (2021) teen-employment application
├── cengiz_run.R                      # Cengiz et al. (2019) minimum-wage application
├── cengiz_application.R              # l0 + DP specification estimators on Cengiz events
├── figures/                          # figures used by the paper
├── *.csv                             # simulation and application outputs
└── archive/                          # earlier prototype code (R/, simulations/)
```

---

## How to run

```r
# Monte Carlo study (variance, coverage, solution paths) -> CSVs + figures
Rscript staggered_sim.R

# Empirical applications
Rscript mpdta_run.R
Rscript cengiz_run.R
```

**Requirements:** R >= 4.0. The applications use the `did` package (Callaway &
Sant'Anna). Build the paper with `latexmk -pdf staggered_did_specification.tex`.

The earlier ℓ₀-only prototype (a smaller DGP and a greedy-only estimator) is kept
in `archive/` for reference.

---

## References

- Callaway, B. & Sant'Anna, P. H. C. (2021). Difference-in-differences with multiple time periods. *Journal of Econometrics*, 225, 200–230.
- Cengiz, D., Dube, A., Lindner, A. & Zipperer, B. (2019). The effect of minimum wages on low-wage jobs. *Quarterly Journal of Economics*, 134, 1405–1454.
- Sun, L. & Abraham, S. (2021). Estimating dynamic treatment effects in event studies with heterogeneous treatment effects. *Journal of Econometrics*, 225, 175–199.
- Wooldridge, J. M. (2025). Two-way fixed effects, the two-way Mundlak regression, and difference-in-differences estimators. *Empirical Economics*, 69, 2545–2587.
