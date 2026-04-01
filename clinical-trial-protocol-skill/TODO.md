# Clinical Trial Protocol Skill — Statistical Calculators TODO

Tracking implementation progress for statistical methods described in
`scripts/AdvancedClinicalTrialProtocolStatistics.md`.

Each calculator ships as a Python script, an R script, and a parity test.

## Tier 1 — Most commonly needed, high regulatory impact

- [x] **Sample size calculator** (t-test, z-test, superiority/non-inferiority)
  - `scripts/sample_size_calculator.py` / `.R`
  - `scripts/tests/test_sample_size_parity.py` — 13 tests
- [x] **Power curve generator** (sweep effect size or sample size)
  - `scripts/power_curve_generator.py` / `.R`
  - `scripts/tests/test_power_curve_parity.py` — 11 tests
- [x] **Sensitivity tables** (2D grid of sample sizes across assumption combos)
  - `scripts/sensitivity_table_generator.py` / `.R`
  - `scripts/tests/test_sensitivity_table_parity.py` — 6 tests
- [ ] **Time-to-event sample sizing** (Schoenfeld + Lachin-Foulkes with accrual modeling)
  - Deps: lifelines (Python), survival (R)
  - Key: power driven by number of events, not patients; piecewise accrual/dropout

## Tier 2 — Frequently needed for pivotal trials

- [x] **TOST equivalence sizing** (Two One-Sided Tests for equivalence margins)
  - `scripts/equivalence_sample_size.py` / `.R`
  - `scripts/tests/test_equivalence_parity.py` — 12 tests
- [ ] **Group sequential interim analysis** (Lan-DeMets alpha-spending, O'Brien-Fleming)
  - Deps: rpact via rpy2 (Python), rpact (R)
  - Key: alpha-spending function, futility boundaries, information fractions
- [ ] **Multiplicity adjustment** (graphical procedures of Bretz et al., Holm, fallback)
  - Deps: statsmodels `multipletests` covers p-value corrections; graphical needs custom impl
  - Key: directed weighted graph of hypotheses with alpha redistribution

## Tier 3 — Advanced, growing regulatory importance

- [ ] **Bayesian assurance** (expected power integrating over effect-size prior)
  - Deps: PyMC or custom Monte Carlo (Python), base R
  - Key: unconditional probability of trial success; always lower than point-estimate power
- [ ] **Adaptive design operating characteristics** (simulation-based)
  - Deps: numpy/scipy + rpact for boundaries
  - Key: sample size re-estimation, combination tests, 10k+ simulated trials
- [ ] **Cluster randomized trial sizing** (design effect DEFF = 1 + (m-1)*ICC)
  - Deps: scipy (Python), base R
  - Key: even small ICCs (0.01-0.05) can nearly double sample size

## Infrastructure

- [ ] **rpy2 bridge to rpact** — single most impactful architectural addition; unlocks group sequential, adaptive, and survival sizing from Python
- [ ] **Missing data / estimand planning** — tipping-point analyses, pattern-mixture models, ICH E9(R1) framework
- [ ] **Stratified randomization planning** — block sizes, Pocock-Simon minimization
