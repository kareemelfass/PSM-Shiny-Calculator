# fct_spec.R : the model specification and its defaults ------------------------
# The app turns its inputs into one `spec` list. Tests and scripts can build the
# same list directly, so every app result can be reproduced without Shiny.
# Default numbers are ILLUSTRATIVE ONLY.

default_econ <- function(arm = c("int", "comp")) {
  arm <- match.arg(arm)
  if (arm == "int") {
    list(upfront = 0, drug_cost = 5000, admin_cost = 300, interval_days = 21, rdi = 0.95,
         c_pf = 2000, c_pd = 4000, c_subseq = 15000, prog_share = 0.8,
         c_terminal = 8000, ae_cost = 1500, ae_qaly = 0.01)
  } else {
    list(upfront = 0, drug_cost = 1500, admin_cost = 300, interval_days = 21, rdi = 1,
         c_pf = 2000, c_pd = 4000, c_subseq = 20000, prog_share = 0.8,
         c_terminal = 8000, ae_cost = 2500, ae_qaly = 0.02)
  }
}

default_spec <- function(arms = list(comp = "Standard care", int = "New drug"),
                         lifetable = NULL, fu_years = NA) {
  hr0 <- hr_to_log(0.7, 0.55, 0.89)
  ep <- function() list(approach = "indep", comp_model = "weibull", int_model = "weibull",
                        hybrid = FALSE, cut = 1, hr = hr0)
  list(
    arms = arms,
    settings = list(horizon = 30, cycle = 7 / 365.25, dr_cost = 0.035, dr_out = 0.035,
                    hcc = TRUE, age = 60, male = 0.5, wtp = 30000, currency = "CU",
                    floor_on = !is.null(lifetable), lifetable = lifetable,
                    pfs_le_os = TRUE, fu_years = fu_years),
    endpoints = list(OS = ep(), PFS = ep()),
    waning = list(type = "none", anchor = "start", start = 5, end = 10, target_hr = 1,
                  frailty = 0, apply_os = TRUE, apply_pfs = TRUE, rationale = ""),
    tot = list(
      comp = list(method = "pfs", model = "weibull", cap = NA, beyond_prog = FALSE),
      int = list(method = "pfs", model = "weibull", cap = 2, beyond_prog = TRUE)
    ),
    econ = list(comp = default_econ("comp"), int = default_econ("int"),
                utility = list(pf = 0.75, pd = 0.60, pf_se = 0.03, pd_se = 0.04)),
    crossover = list(comp_os_basis = "itt", subseq_includes_int = TRUE)
  )
}

# Read IPD csv and convert time to years
read_ipd <- function(path, unit = "months") {
  d <- utils::read.csv(path, stringsAsFactors = FALSE)
  names(d) <- tolower(names(d))
  need <- c("arm", "endpoint", "time", "event")
  miss <- setdiff(need, names(d))
  if (length(miss)) stop("Missing column(s): ", paste(miss, collapse = ", "))
  d$endpoint <- toupper(trimws(d$endpoint))
  d$endpoint[d$endpoint == "TOT"] <- "ToT"
  bad <- setdiff(unique(d$endpoint), c("OS", "PFS", "ToT"))
  if (length(bad)) stop("Unknown endpoint value(s): ", paste(bad, collapse = ", "),
                        ". Use OS, PFS or ToT.")
  if (!all(d$event %in% c(0, 1))) stop("event must be 0 (censored) or 1 (event).")
  d$time <- as.numeric(d$time) * time_factor(unit)
  d$arm <- as.character(d$arm)
  # Optional patient id: links a patient's OS, PFS and ToT rows (paired bootstrap)
  d$id <- if ("id" %in% names(d)) {
    x <- trimws(as.character(d$id))
    x[x == ""] <- NA
    ifelse(is.na(x), NA_character_, paste(d$arm, x, sep = "::"))
  } else NA_character_
  d[d$time > 0 | d$event == 0, c(need, "id")]
}

# TRUE when every row has a patient id (bootstrap keeps OS, PFS and ToT together)
ipd_is_paired <- function(ipd) "id" %in% names(ipd) && nrow(ipd) > 0 && all(!is.na(ipd$id))

read_lifetable <- function(path) {
  lt <- utils::read.csv(path)
  names(lt) <- tolower(names(lt))
  if (!all(c("age", "qx_male", "qx_female") %in% names(lt))) {
    stop("Life table needs columns: age, qx_male, qx_female (annual probability of death).")
  }
  lt[order(lt$age), c("age", "qx_male", "qx_female")]
}

# Fit everything a spec needs (used by tests and scripts)
fit_all <- function(ipd, arms, engine = fit_engine()) {
  list(
    OS = fit_endpoint(ipd, "OS", arms, engine),
    PFS = fit_endpoint(ipd, "PFS", arms, engine),
    ToT = if (any(ipd$endpoint == "ToT")) fit_endpoint(ipd, "ToT", arms, engine) else NULL
  )
}

# Script entry point: spec + fits -> everything
run_model <- function(spec, fits) {
  t <- make_grid(spec$settings$horizon, spec$settings$cycle)
  pieces <- predict_pieces(spec, fits, t)
  curves <- assemble_curves(pieces, spec)
  res <- run_engine(curves, spec)
  list(pieces = pieces, curves = curves, res = res,
       validation = validate_model(spec, pieces, curves, res, fits))
}
