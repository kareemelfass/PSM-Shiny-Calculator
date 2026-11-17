# Partitioned survival model (Shiny), V3

## What is new in V3: treatment-effect waning (Doc 12)

- **Four structures, labelled with their appraisal role:** sustained effect (needs evidence of durable
  benefit), equal hazards from a time point (pragmatic base-case candidate; survival already gained is
  kept), gradual equalisation (HR moves linearly to a target), stop and drop (extreme stress test).
- **Timing anchor:** years since treatment start, or years after the maximum treatment duration
  (NICE TA650: effect ceasing 3 years after a 2-year treatment stops = year 5).
- **Target HR after waning** (default 1) for partial-waning scenarios.
- **Individual-level waning (frailty):** with a gamma frailty variance > 0, the conditional
  (individual) effect wanes and the population-level HR goes above the target, because surviving
  intervention patients are frailer. Constraining the marginal HR to 1 overestimates benefit and
  underestimates the ICER (Jennings et al. 2024). Default 0 = V2 behaviour.
- **Waning explorer:** every structure x a set of start times in one run, with incremental QALYs,
  ICER, landmark OS and whether the conclusion changes (NICE TA798, TA1099: reviewers test timings and
  convergence structures). Warns when results are stale; exported with a `matches_current_inputs` flag.
- **Plots:** timeline (treatment, evidence, effect), intervention OS under each structure, absolute
  OS advantage over time (equal hazards are not equal survival), implied HR and hazards.
- **DSA:** waning start (and duration, if gradual) in the tornado. **PSA:** optional uniform sampling
  of the waning start.
- **Validation:** sustained effect after a stopping rule; waning starting inside observed follow-up;
  stop-and-drop instant loss (size reported: a survival jump, not a finite hazard); anchor without a
  treatment cap; post-waning marginal HR under frailty; waning that lowers the intervention's hazard
  (weighted by patients at risk); missing justification text.
- Defaults reproduce V2 exactly (all V2 tests pass unchanged). Tests reproduce Doc 12's worked example:
  76.3% vs 58.3% at year 3; 40.7%, 21.7% and 16.5% at year 10; gap 18.1, 8.8, 5.1 points.

### Waning maths

Per cycle, on the chosen scale (population hazards, or individual hazards if frailty > 0):

    h_int,waned(t) = (1 - w(t)) x h_int(t) + w(t) x target_HR x h_comp(t)

w(t) = 0 before the start, rises linearly to 1 by the end (gradual) or jumps to 1 (equal hazards).
Frailty: conditional cumulative hazard Hc = (exp(theta Hm) - 1) / theta, blended, then mapped back
with Hm = log(1 + theta Hc) / theta. Stop and drop copies the comparator curve from the start time.

References (from Doc 12): Taylor et al. PharmacoEconomics 2024;42:1181-1196. Jennings et al. Value
Health 2024;27:347-355. Latimer, NICE DSU TSD 14 (2013). NICE TA650, TA798, TA1099.

# V2 notes

## What is new in V2

- **Survival-curve uncertainty in the PSA** (PSA tab, "Survival curve uncertainty"):
  - *Parameter covariance*: each fitted model's parameters drawn from N(estimate, covariance).
    survreg: coefficients and log(scale) with `vcov()`; flexsurv: `normboot.flexsurvreg()`.
    Both arms of a joint model use the same draw. Fast.
  - *Bootstrap*: patients resampled within arm and the chosen models refitted each run.
    With an `id` column, a patient's OS, PFS and ToT rows are resampled together, which keeps
    their correlation (Doc 13 Case 16). Slower: use 500 runs and parallel workers.
  - *None*: curves fixed (V1 behaviour; reproduces V1 draws exactly for the same seed).
- Hybrid KM + model curves: the KM part is re-estimated in the bootstrap, fixed under parameter sampling.
- Fitted ToT curves are sampled too; HRs entered by the user are sampled as in V1.
- Parallel workers (base R `parallel`, works on Windows). Results are identical to a serial run.
- New PSA outputs: survival uncertainty bands vs the deterministic curves, deterministic vs PSA mean
  table, convergence plot, and a table of invalid runs with reasons (excluded, never fixed).
- New validation rows and Excel sheets: PSA_settings, PSA_survival_bands, PSA_vs_deterministic.
- Optional `id` column in the patient data; the synthetic example now has one.
- flexsurv sampling has a built-in self-check: it stops (suggesting the bootstrap) if the
  parameter route does not reproduce the fitted curve.


Three-state partitioned survival model (progression-free, progressed, dead) for HTA work.
All maths lives in plain R functions in `R/fct_*.R`; the Shiny modules (`R/mod_*.R`) only collect
inputs and display results. Every number the app shows can be reproduced from a script.

**All example data, the life table and default costs/utilities are synthetic or illustrative.**

## Run

```r
install.packages(c("shiny", "bslib", "survival", "flexsurv", "ggplot2", "writexl", "testthat"))
shiny::runApp("C:/path/to/PSM_Shiny_V3")   # the folder that contains app.R
```

Or open `app.R` in RStudio and click **Run App**. `app.R` loads its own `R/` files, so it also works
when sourced or run line by line, as long as it can find its folder (it checks the sourced file's
location, the working directory, then the file open in RStudio).

On Windows, "Extract all" can create `PSM_Shiny_V3/PSM_Shiny_V3/`: point `runApp()` at the folder that directly
contains `app.R`.

Requires R >= 4.1, shiny >= 1.7, bslib >= 0.6. Tested here with R 4.3.3, shiny 1.8.0, bslib 0.6.1,
flexsurv 2.3.3 and survival 3.5-8. No new packages in V2 (`parallel` ships with R).
For submission work: `renv::init()` then `renv::snapshot()` in the app folder.

Without flexsurv the app still runs, using `survival::survreg` (exponential, Weibull, log-normal,
log-logistic). With flexsurv it adds Gompertz, gamma, generalised gamma and Royston-Parmar splines (1 to 3 knots).

## First thing to do on your machine

```r
testthat::test_dir("PSM_Shiny_V3/tests/testthat")
```

All tests passed in the build environment with flexsurv installed, including: flexsurv and survreg
give the same fits and the same Weibull PSA uncertainty; every flexsurv model's sampled curves are
centred on its fitted curve; bootstrap keeps patients whole; parallel equals serial.

## Script use (no Shiny)

```r
for (f in list.files("PSM_Shiny_V3/R", "^fct_|^utils", full.names = TRUE)) source(f)
library(survival)
ipd  <- read_ipd("PSM_Shiny_V3/example/example_ipd_months.csv", "months")
lt   <- read_lifetable("PSM_Shiny_V3/example/illustrative_lifetable.csv")
arms <- list(comp = "Standard care", int = "New drug")
fits <- fit_all(ipd, arms)
spec <- default_spec(arms, lifetable = lt, fu_years = 4)
out  <- run_model(spec, fits)
out$res$icer; out$validation
```

## Data format

CSV with columns `arm, endpoint, time, event` and optional `id`. `endpoint` is OS, PFS or ToT; `event` is 1
(event) or 0 (censored). `id` links one patient's rows across endpoints (used by the PSA bootstrap).
Choose the time unit on upload. The intervention arm may have no data when its curves come from a hazard ratio
(for example an NMA). Digitised curves can be converted to patient data on the Data tab (Guyot et al. 2012).

## Methods and where they come from

| Feature | Implementation | Source document |
|---|---|---|
| State occupancy | PF = PFS; PD = OS - PFS; Dead = 1 - OS | Doc 8; TSD 19 |
| Hazard ratio | S_int(t) = S_comp(t)^HR; HR direction converted to intervention vs comparator; SE from the 95% interval on the log scale | Doc 3; Doc 9 B4 |
| Model fit | AIC, BIC (n = patients), dAIC/dBIC within each dataset; KM plus fits; extrapolated hazards | Doc 6 |
| Hybrid KM + model | KM to the cut-point, then S_model(t) / S_model(cut) x KM(cut); risk set shown | Doc 6 |
| PH diagnostics | Cox HR, Schoenfeld test and residual plot, log-cumulative hazard plot | Doc 3 |
| Plausibility | Landmarks vs expert benchmark ranges; candidate-model table (AIC, landmarks, mean); hazard vs general population | Doc 7 |
| Waning | Hazard blend to a target HR; equal, gradual, stop and drop; timing from start or treatment stop; optional individual-level (frailty) waning; explorer, DSA, PSA | Doc 12; Doc 9 C4a |
| General-population floor | Life-table hazard (mixed by sex) as a lower bound on OS and PFS hazards | Doc 7 |
| Time on treatment | Fitted ToT or equal to PFS; maximum duration; beyond progression on/off; always at or below OS | Doc 11 |
| Drug cost | Area under ToT x cost per administration x (365.25 / interval) x relative dose intensity | Doc 11; Doc 13 |
| Crossover | Comparator OS basis (ITT or adjusted) must match subsequent-therapy costs; mismatch flagged | Doc 10 |
| Validation | Full-horizon monotonicity, 0 to 1, occupancy sums to 1, HR = 1 reproduces reference, totals, constraints reported (never silent) | Doc 13 Cases 14, 16 |
| PSA | log-HR normal (OS/PFS correlation optional), utilities beta, costs gamma; survival curves by parameter covariance or paired bootstrap; per-run validity with reasons; ICER = ratio of means; convergence | Doc 13 Cases 16, 17 |
| Extrapolation dependence | Share of incremental life-years after trial follow-up | Doc 9 C4a |

Discounting: 1 / (1 + r)^t at each cycle mid-point. Half-cycle correction: trapezoid rule.
New progressions per cycle = fall in PFS x share of PFS events that are progressions (an input; set it from the trial).

## Not in this version

- Under parameter sampling, OS, PFS and ToT models are drawn independently (use the paired bootstrap
  to keep their correlation). Joint OS/PFS models (e.g. copulas, multi-state) are not included.
- The PSA runs in the foreground (progress bar); the app is busy until it finishes.
- Crossover adjustment methods (RPSFTM, IPCW, two-stage): supply already-adjusted data instead.
- Cure and mixture-cure models; state transition model cross-check (TSD 19, recommendation 11).
- Age-adjusted utilities.

## Structure

```
app.R                  page layout and wiring
R/fct_survival.R       fitting adapter (flexsurv or survreg), KM, hybrid, Cox
R/fct_effects.R        HR, waning, general-population floor, constraints
R/fct_engine.R         occupancy, life-years, QALYs, costs, ICER
R/fct_validation.R     automated checks, landmarks
R/fct_uncertainty.R    DSA, PSA driver and PSA summaries
R/fct_psa_surv.R       V2: survival samplers (parameter covariance, bootstrap), parallel workers
R/fct_waning.R         V3: waning structures, diagnostics, structure comparison, explorer
R/fct_guyot.R          IPD reconstruction from digitised KM
R/fct_spec.R           defaults, file readers, run_model()
R/mod_*.R, R/plots.R   Shiny modules and plots
tests/testthat/        known-answer tests (engine, floor, Guyot, PSA, V2 survival sampling, V3 waning vs Doc 12)
example/               synthetic IPD, digitised curve, illustrative life table
data-raw/              script that generated the example files
```

## PSA from a script

```r
t      <- make_grid(spec$settings$horizon, spec$settings$cycle)
pieces <- predict_pieces(spec, fits, t)
runs   <- run_psa(pieces, spec, n = 1000, seed = 2026, surv_method = "boot",
                  fits = fits, ipd = ipd, cores = 4, code_dir = "PSM_Shiny_V3/R")
psa_summary(runs, wtp = 30000)
psa_vs_deterministic(runs, run_engine(assemble_curves(pieces, spec), spec))
```
