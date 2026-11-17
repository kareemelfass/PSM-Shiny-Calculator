# make_example.R : builds the synthetic example files in example/ -------------
# ALL DATA ARE SIMULATED FOR TESTING. They describe no real drug or trial.
# Run from the app folder:  Rscript data-raw/make_example.R

set.seed(2026)
n <- 250
accrual <- 1.5          # years of recruitment
study_end <- 4          # years from first patient: maximum follow-up = 4 years

sim_arm <- function(arm, n, os_shape, os_scale, prog_rate, stop_rate, beyond_prog) {
  entry <- stats::runif(n, 0, accrual)
  admin_cens <- study_end - entry
  dropout <- stats::rexp(n, 0.03)
  cens <- pmin(admin_cens, dropout)

  os_true <- stats::rweibull(n, shape = os_shape, scale = os_scale)
  prog_true <- stats::rexp(n, prog_rate)
  pfs_true <- pmin(os_true, prog_true)
  stop_tox <- stats::rexp(n, stop_rate)
  extra <- ifelse(stats::runif(n) < beyond_prog, stats::rexp(n, 3), 0)
  tot_true <- pmin(os_true, pmin(pfs_true + extra, stop_tox))

  mk <- function(endpoint, true) {
    data.frame(id = sprintf("%s%03d", substr(gsub(" ", "", arm), 1, 3), seq_len(n)),
               arm = arm, endpoint = endpoint,
               time = round(pmin(true, cens) * 12, 3),      # months
               event = as.integer(true <= cens))
  }
  rbind(mk("OS", os_true), mk("PFS", pfs_true), mk("ToT", tot_true))
}

ipd <- rbind(
  sim_arm("Standard care", n, os_shape = 1.15, os_scale = 1.9,
          prog_rate = 0.95, stop_rate = 0.25, beyond_prog = 0.05),
  sim_arm("New drug", n, os_shape = 1.05, os_scale = 2.9,
          prog_rate = 0.60, stop_rate = 0.30, beyond_prog = 0.25)
)
dir.create("example", showWarnings = FALSE)
utils::write.csv(ipd, "example/example_ipd_months.csv", row.names = FALSE)

# Digitised-curve example (comparator OS) for testing Guyot reconstruction
os_c <- ipd[ipd$arm == "Standard care" & ipd$endpoint == "OS", ]
km <- survival::survfit(survival::Surv(time, event) ~ 1, data = os_c)
digi <- data.frame(time = c(0, km$time), surv = round(100 * c(1, km$surv), 2))
digi <- digi[!duplicated(digi$surv), ]
risk_t <- seq(0, 42, by = 6)
risk <- data.frame(time = risk_t,
                   n_risk = summary(km, times = risk_t, extend = TRUE)$n.risk)
utils::write.csv(digi, "example/example_digitised_os_standard_care.csv", row.names = FALSE)
utils::write.csv(risk, "example/example_nrisk_os_standard_care.csv", row.names = FALSE)

# Illustrative life table (Gompertz-Makeham shape). NOT real population data:
# replace with a national life table (e.g. WHO Global Health Observatory) for any analysis.
age <- 0:100
lt <- data.frame(
  age = age,
  qx_male = round(1 - exp(-(0.0008 + 0.000040 * exp(0.093 * age))), 6),
  qx_female = round(1 - exp(-(0.0006 + 0.000022 * exp(0.096 * age))), 6)
)
utils::write.csv(lt, "example/illustrative_lifetable.csv", row.names = FALSE)
cat("Example files written to example/\n")
