# PAmeasures 2.0.0

## New features

* Extend performance assessment to right-censored survival outcomes, competing
  risks, case-cohort designs and nested case-control designs. Prediction and
  evaluation are separate, allowing externally fitted predictions to be used.
* Provide twenty primary `pa_` interfaces for prediction, evaluation,
  calibration, model summaries, sampling weights, plotting and simulation.
  Survival methods are adapted from TimeMetric; its historical release notes
  are retained below.
* Add `pa_nlm()` for continuous-outcome R-squared and L-squared, returning
  unrounded numeric results. Retain `pam.nlm()` with its historical field
  names and four-decimal character format, using the same calculation core.
* Provide integrated and point-in-time R-squared measures, concordance,
  Brier scores and time-dependent AUC for their supported outcome settings.
* Provide ICI, E50 and E90 calibration assessment at a specified time for
  right-censored and competing-risk cohorts, and sampling-weighted calibration
  for right-censored two-phase designs. Evaluation defaults include ICI;
  E50 and E90 can be requested explicitly.

## Compatibility and input handling

* All primary function names use `pa_`. `pam.nlm()` is the sole supported
  historical public function name. The PAmeasures 0.1.0 names `pam.censor()`,
  `pam.coxph()` and `pam.survreg()`, the `tm_` interfaces and earlier TimeMetric
  script aliases are not exported. Full API compatibility with 0.1.0 is not
  claimed. The `moore` dataset is retained.
* Validate paired finite vectors explicitly for continuous outcomes. Constant
  predictions use an intercept-only correction. Undefined denominators return
  NA with a warning while retaining any other defined measure. Inputs causing
  overflow or underflow in squared deviations require rescaling. These
  documented edge rules replace historical incidental NaN/Inf behavior.
* Retain the GPL-3 license, original PAmeasures authorship and TimeMetric MIT
  notices for the adapted source.

## Documentation and analysis workflows

* Update package help, examples and citation metadata for the unified API.
* Include repository scripts for the paper simulations, real-data analyses
  and manuscript figures. A unified runner checks dependencies and runs all
  analyses, with simulation checkpoints for resumption. Manuscript refreshes
  can redraw saved results without refitting models. Generated results and
  clinical predictions are excluded from the distributed source package.
* Require a new simulation output directory when computation sources,
  settings or package versions change, preventing incompatible checkpoints
  from being mixed.

# TimeMetric 0.2.1

* Correct cause-specific AFT CIF integration near time zero using adaptive quadrature on the target-event probability scale. Fixed-scale Weibull models with scale greater than one no longer lose the initial probability mass.
* Align the Framingham analysis scripts with the manuscript's 12-year horizon and correct the CIF time-summary documentation.

## Calibration

* Three- and four-valued ordinary risks now use linear Cox calibration,
  matching unit-weighted calibration without spline errors or console output.
* Saved-result upgrades defer manuscript publication until PBC, Framingham
  and GBSG all contain ICI. Explicit input paths work without a published-run
  manifest. A simulated end-to-end migration test preserves prior metrics,
  original files and the active manifest until the upgrade is complete.

* Added `tm_calibration()` for right-censored event probabilities and
  `tm_calibration_cr()` for a target-cause CIF, both at one specified horizon
  and without requiring the original prediction model. Results contain ICI,
  E50 and E90 plus calibrated probabilities in original subject order.
* Right-censored calibration uses `pmcalibration` restricted cubic splines;
  competing-risks calibration uses `prodlim` nearest-neighborhood
  Aalen--Johansen estimates. Both use `pmcalibration::cal_metrics()` for the
  summaries, retaining repeated predictions with their subject frequencies.
  `pmcalibration` and `prodlim` are now declared runtime dependencies.
* Ordinary and competing-risks evaluators and comparison tables include ICI
  by default. E50/E90 are opt-in (`metrics = c("ici", "e50", "e90")` or
  `"all"`). `tm_fit_and_eval()` now defaults to `metrics = NULL` for this
  same ICI-only calibration selection; explicit `"all"` includes all three.
  `calibration_args` forwards knot/bandwidth settings to the new functions.
* Calibration uses the requested outcome horizon and the last supplied
  prediction-grid value at or before it. Direct probability interfaces and
  fit-and-evaluate can use predictions calculated exactly at the horizon.
  Zero/one probabilities are never silently dropped; ordinary survival
  fitting clips them inward with a warning and summaries use the originals.
  Non-estimable calibration fits warn and leave other metrics available.
  Calibration confidence intervals are not included.

* Right-censored two-phase evaluation now includes ICI by default, with
  E50/E90 available explicitly or via `"all"`. `tm_calibration(case_weights = ...)`
  fits a sampling-weighted Cox spline model and computes weighted summaries;
  `tm_evaluate_two_phase()` and `tm_sample_design()` forward aligned sampling
  weights and `calibration_args`. Calibration-only calls need no original
  prediction model, mean-time predictions, or separate censoring fit.
  Weighted constant predictions use Kaplan-Meier; all fits retain original
  boundary probabilities for metrics. `Hmisc` is a direct dependency for the
  shared spline basis and weighted quantiles. Design-based intervals and
  sampling-adjusted competing-risks calibration are not implemented.

## Paper integration

* Simulation 1 now uses common random numbers: all 16 scenarios share the
  training seed 1234 (as in Simulation 2) and test seed r for repetition r,
  replacing scenario-specific seeds. Its Cox training fits use
  `timefix = FALSE`, because with 1/v = 2 the default merged thousands of
  distinct times into ties and shrank the beta = 1.5 coefficients by about
  16%. Simulation 3 keeps the Simulation 1 beta = 1, v = 1 model, which is now
  the seed-1234 fit.

* Figure 6 is drawn as a final figure when its only missing values are
  Simulation 3 ICI repetitions with no follow-up at the horizon. Those boxes
  use the finite repetitions; the cells and counts are written to
  `paper.fig6.partial_repetitions.csv` (n = 300, 70% censoring, 90% horizon:
  52/100) for the manuscript caption.

* Simulation ICI now requests fitted-model probabilities at exactly its outcome
  horizon, rather than the preceding test-follow-up grid point. Other metrics
  retain their existing evaluation procedures. The Simulation 3 reference
  remains the mean across the uncensored n=3000 repetitions. A paired replay
  script compares saved and corrected ICI without refitting prediction models.

* Simulation 1, both Simulation 2 parts and Simulation 3 now include ICI at
  their existing evaluation horizons, with five-knot Cox spline calibration.
  ICI is stored without rounding and enters replication summaries, Simulation
  3 references/errors, figure panels and report tables. Calibration diagnostics
  are saved per job and model. Resume and report checks require ICI rows.
  Existing simulation outputs must be kept in a separate directory.

* Competing-risks C-index now uses reciprocal censoring-survival products,
  following Wolbers et al. (2014), equation (3.4). Previously it multiplied
  comparable pairs by those products. Only contributing pairs are inverted,
  avoiding undefined weights in unused rows after the evaluation horizon.
  Existing competing-risks C-index results require reevaluation.
* Point R-squared now includes target events at the evaluation time (`T <= t`)
  in both survival and competing-risks evaluation, matching the manuscript.
  Ordinary survival point R-squared compares the event indicator with `1 - S`
  rather than `S`; this changes its adequacy and pseudo-R-squared values.
* Competing-risks generalized R-squared now uses the event-observation IPCW
  strategy in `PA measure new/pam.survreg.update.R`: original follow-up times
  estimate censoring and evaluate weights, and only observed events of any
  cause contribute. Both interval and point versions expose Cox (default,
  optionally conditional) and KM censoring estimates, with a probability floor
  of 1e-4. Late censoring remains zero-weight by design. The adapted Cox
  baseline/predictor centering and tied-time left limits are kept consistent.
  These changes alter CR R-squared results; saved tables need recomputation.
* Cox predictions now obtain the baseline hazard exclusively from the training
  fit. Validation labels cannot change survival probabilities on a fixed grid.
  Training weights, centering, ties, factors and strata follow `survival`.
  The legacy evaluation-grid integration remains explicit; this repair does
  not claim that unrestricted mean survival is identified.
  Offset models use the same weighted training-offset centering as the baseline;
  reference tests cover in-sample and validation predictions, case weights,
  strata, and Efron, Breslow and exact tie methods.
* Matched NCC weights exclude cases from other matching strata. Each case is
  a sampling occasion, including tied cases; small risk sets provide all
  available controls. `m` remains the number of controls per case, with
  full-cohort input. Simultaneous cases are excluded from control draws and
  events precede censoring at ties. Paper examples and tests use this convention.
* IPCW weights for ordinary survival R-squared, point R-squared and the
  Cox/AFT prediction helpers now use one censoring left limit G(t-) for all
  rows tied at t. Previously, results with tied times depended on row order,
  and point R-squared was off when a subject was censored exactly at t*.
  Data with tied times can give slightly different values.
* The competing-risks Cox censoring model uses Breslow ties, so the
  no-covariate model equals exp(-Nelson-Aalen) with tied censoring times.
* `tm_predict_coxph()` supports strata-by-covariate interactions such as
  `x * strata(g)` when `new_data` is given, and no longer repeats survival's
  interaction warning on every call.
* `tm_fit_and_eval()` computes R-squared only when requested. Cox/AFT
  R-squared with undefined censoring weights returns missing values with a
  suppressible diagnostic instead of recycling shorter calibration vectors.
* Added independent prediction/inclusion-probability checks and paper analysis
  scripts. PBC is compared before/after within 0.2.1; full simulation execution
  is an explicit separate action. See the paper workflow in `README.md`.

## CRAN review corrections

Resubmission addressing the CRAN review of 0.2.0 and correcting numerical
errors found during independent validation. Several metric values change;
analyses using previous versions should be rerun.

* `DESCRIPTION`: the Description no longer starts with the package name,
  quotes software names, spells out AUC, and cites the methods it implements
  with DOIs. The predecessor package is now correctly named 'PAmeasures'.
* `tm_predict_coxph()` and `tm_predict_survreg()` default to `predict = TRUE`
  rather than `predict = T`; `T`/`F` are no longer used anywhere in the code.
* Diagnostics are no longer written with `print()` or `cat()`. A failed model
  fit inside `tm_predict_cif()` and a failed weighted least-squares fit inside
  `tm_predict_survreg()` now raise a warning.
* The time-dependent AUC in `tm_fit_and_eval()` now calls
  `survivalROC::survivalROC()` directly instead of shipping a copy of that
  function's GPL-licensed source, and 'survivalROC' is added to Imports.
* `tm_simulate_fine_gray()` and `tm_plot_summary()` no longer reset the random
  number generator by default. Their `seed` argument now defaults to `NULL`;
  pass a seed to get the previous reproducible behavior.
* 'tidyr' is required in Imports because the default two-phase AUC uses it
  through 'yardstick'. A CI job tests the default metrics using only required
  dependencies. The unused 'expint' dependency has been removed.
* `tm_fit_and_eval()` excludes incomplete training/evaluation rows with a
  suppressible warning before invoking any metric backend. This also prevents
  unsuppressible missing-value output from 'survivalROC'.

## Numerical and API corrections

* Competing-risk `td_auc` can now be requested without also requesting
  `brier_score`.
* `tm_fit_and_eval()` matches values to metric names when a subset is supplied
  in a different order. Previously it could silently label a Brier score as
  concordance (or vice versa). Unsupported models and metrics now give clear
  errors.
* `tm_predict_survreg()` retains the covariate column when predicting from a
  data frame containing a single predictor.
* Cox survival probabilities now come from the fitted training model through
  `survival::survfit()`. Previously the baseline was re-estimated using
  validation outcomes, leaking information and ignoring fitted case weights.
* Finite-horizon Weibull and lognormal AFT means now compute the integral of
  survival, including the probability mass beyond the restriction time.
  Prediction objects retain original follow-up times and event labels.
* Explained-variation calculations reorder predictions and sampling weights
  together with observed times.
* Point pseudo R-squared now compares event indicators with event probabilities
  (`1 - survival`), including events at the evaluation horizon. Point
  calibration-component and pseudo R-squared values consequently change.
* `tm_fit_and_eval()` orients AFT scores by their model meaning and retains
  AUC values below 0.5. Previously `max(AUC, 1 - AUC)` could make reversed
  held-out predictions appear highly discriminative.
* NCC weights interpret `m` as the documented number of controls per case.
  Sampling opportunities are restricted to cases in the matching stratum;
  tied failures and censoring follow the documented risk-set convention.
  NCC metric values and finite-horizon AFT results change as a consequence.
  Previously reported manuscript results need to be regenerated.

# TimeMetric 0.2.0

Package-quality release preparing resubmission to JOSS. The evaluation
methodology is unchanged; metric *values* are identical to 0.1.0 except where a
defect is listed below.

## Installation

* **The package can now be installed.** 0.1.0 imported `randomForestSRC` in
  `NAMESPACE` without declaring it anywhere in `DESCRIPTION`, so
  `R CMD INSTALL` failed on any machine that did not already have it and
  `remotes::install_github()` could not succeed.
* All runtime dependencies are declared: `rms`, `survminer`, `expint`, `pec`,
  `magrittr`, `purrr`, `tibble` and `utils` were used but undeclared.
* `randomForestSRC` and `cmprsk` moved to `Suggests` as optional backends of
  `tm_predict_cif()`, each guarded with an informative error naming the package
  to install.
* **`rms` is no longer a dependency.** It was used only to fit a model for the
  since-withdrawn `r_sh` estimator and to extract linear predictors, both of
  which are now `survival` calls. Because `rms` requires R >= 4.4.0, removing
  it lowers the declared minimum from R 4.4.0 to **R >= 4.1.0**, the floor set
  by `survival` and `dplyr`.

## Renamed API

Every exported function now carries a `tm_` prefix, correcting the `survial` and
`surverg` misspellings and dropping the `pam.` prefix inherited from PAmeasures.
All previous names still work and forward to their replacement with a
deprecation warning; see `?"TimeMetric-deprecated"`.

* `tm_evaluate_two_phase()` is **newly exported**. The case-cohort and nested
  case-control functionality was previously internal and unreachable by users.
* `tm_metric_names()` returns the canonical metric identifiers.

## Removed metrics

* Removed `r_sh` and `r_e` from the pre-release API following a maintainer
  decision to exclude metrics that are not sufficiently validated for the first
  CRAN release. Requesting either name, or any of the legacy spellings `R_sh`,
  `R_E` and `R_sph`, now raises an error naming the withdrawal.
* `tm_metric_names()` consequently returns 11 canonical names rather than 13.
* Every other metric is numerically unchanged, verified against values recorded
  before the removal across all four public evaluators, including the
  competing-risks and two-phase paths.

## Bug fixes

* **Point-in-time metrics could differ between platforms, and their values have
  changed.** The evaluation time was chosen with
  `which.min(abs(event_time - t_star))`. When `t_star` is the median of an
  even-length vector the two middle observed times are mathematically
  equidistant from it, so the winner was decided by a rounding
  difference in the last bit of the result -- which falls differently on platforms with 80-bit `long double`
  (x86_64 Linux and Windows) than on those where it is plain double (arm64
  macOS). The same data therefore produced different `pseudo_r2_point`,
  `r2_point`, `l2_point`, `brier_score` and `td_auc` on different machines.

  The rule is now nearest observed time with ties broken toward the **earlier**
  event time, which is also independent of input row order. All metrics now
  agree bit-for-bit across arm64 macOS, x86_64 Linux and x86_64 Windows.
  Integrated and concordance metrics were never affected.

* **The Brier score's censoring weights were wrong, and `brier_score` from
  `tm_fit_and_eval()` changes on data with tied times.** `Gt()`, the
  Kaplan-Meier estimate of the censoring distribution behind the
  inverse-probability-of-censoring weights, had three independent defects: it
  derived interpolation indices from a sorted table but applied them to the raw
  unsorted input, so G(t) was neither monotone nor invariant to row order; it
  identified censoring events as `status == min(status)`, which inverts on data
  with no censoring and made G decay from 1 instead of staying at 1; and it
  silently replaced a zero censoring survival with the smallest positive value,
  concealing an undefined weight.

  `Gt()` is now a reverse-Kaplan-Meier **step function** with no interpolation
  between jump times, using G(t) for the evaluation-time weight and G(t-) for
  subject-specific event-time weights, and the reverse-Kaplan-Meier tie
  convention. It returns `NA` beyond the last observed time and **errors** when
  the censoring distribution is exhausted, rather than substituting a value.
  Values agree exactly with `pec::ipcw()`.

  Reported Brier scores on tied data at an evaluation time between jump times
  were substantially overstated: on one 200-observation fixture the reported
  value was 0.75 against a correct 0.22. Metrics other than `brier_score` are
  unaffected, as is `tm_survival_eval()`, which computes the Brier score by a
  different route.

* `tm_fit_and_eval()` (was `pam.survival_eval()`) **never worked**. It passed
  `covariates=`/`newdata=` to functions taking `covs=`/`new_data=`, so every
  call failed with `unused arguments`.
* `tm_plot_pred()` (was `plot_pred()`) returned a **silently empty plot** on its
  documented default call. `sample_index` defaulted to `NULL` and the function
  then subset with `x[NULL]`. `NULL` now means "plot every subject".
* `survival::concordancefit()` was called without being imported, so the primary
  evaluation function failed under `library(TimeMetric)` alone. It worked only
  if the user happened to have attached `survival` separately.
* **`tm_survival_eval_cr()` rejected every legacy metric spelling.** It carried
  a duplicate metric-validation block that ran *before* name normalization, so
  `"Harrells_C"`, `"C_index"` and `"Brier Score"` all failed on the
  competing-risks path with `Invalid metrics:` while working on every other
  entry point. Validation now runs after normalization, as it does elsewhere,
  and the documented legacy spellings are accepted there too. Genuinely unknown
  names are still rejected.
* `tm_survival_eval_cr()` returned its `Value` column as **character** while
  `tm_survival_eval()` returned numeric, so `is.finite()` was silently `FALSE`
  for every competing-risks value. Both now return numeric.

## Metric names

Standardized to one ASCII identifier per quantity. Legacy spellings are still
accepted, case-insensitively and tolerant of apostrophe style, with a
deprecation warning:

`brier_score`, `c_index`, `harrell_c`, `l_square`, `l2_point`, `pseudo_r2`,
`pseudo_r2_point`, `r_square`, `r2_point`, `td_auc`, `uno_c`

* `Harrell's C` and `Uno's C` previously required a U+2019 curly apostrophe to
  match, which made them close to impossible to discover.
* The pseudo R-squared family carried four spellings, one of them the misspelled
  `Pesudo_R`. `pseudo_r2` and `pseudo_r2_point` remain distinct measures.
* Time-dependent AUC carried different labels in the two-phase evaluator and its
  own summary wrapper, so results could not be joined. Both emit `td_auc`.

## Removed

Thirteen unreachable functions were deleted after a call-graph analysis and, for
the metric implementations, an equivalence audit against the original
Stare-Perme-Henderson reference code:

* `pam.rsph_metric()` and `pam.rsh_metric()` were duplicate implementations of
  the two metrics that have since been withdrawn entirely (see **Removed
  metrics**).
* `pam.Brier_metric()` duplicated the Brier score already computed by a
  weight-aware implementation.
* `pam.prediction_survial_eval()`, `pam.prediction_metrics()`,
  `pam.prediction_metrics_cr()`, `pam.coxph()`, `pam.nlm()`, `pam.survreg()`
  and related helpers were unreachable from any export.

## Testing and infrastructure

* New `testthat` suite: 498 assertions across 100 tests, none skipped.
* Every surviving metric is pinned to values recorded before the `r_sh` / `r_e`
  withdrawal, across all four public evaluators including the competing-risks
  and two-phase paths.
* GitHub Actions: `R-CMD-check` across five platform/version combinations, a
  full-dependency job that fails if any test skips, test coverage, and the
  existing JOSS draft-PDF build, whose path trigger had pointed at a filename
  that never existed.
* `R CMD check --as-cran` reports 0 errors and 0 warnings.
