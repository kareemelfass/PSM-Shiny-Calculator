# fct_uncertainty.R : one-way DSA and PSA ----------------------------------------
# Survival model parameters are fixed in this version. PSA samples log-HRs
# (where the HR approach is used), utilities and costs.

# Methods-of-moments helpers
beta_pars <- function(m, se) {
  k <- m * (1 - m) / se^2 - 1
  c(a = m * k, b = (1 - m) * k)
}
rbeta_ms <- function(n, m, se) {
  if (!is.finite(se) || se <= 0 || m <= 0 || m >= 1) return(rep(m, n))
  p <- beta_pars(m, se)
  if (any(p <= 0)) return(rep(m, n))
  stats::rbeta(n, p[["a"]], p[["b"]])
}
rgamma_mcv <- function(n, m, cv) {
  if (!is.finite(m) || m <= 0 || cv <= 0) return(rep(m, n))
  stats::rgamma(n, shape = 1 / cv^2, scale = m * cv^2)
}

cost_fields <- c("upfront", "drug_cost", "admin_cost", "c_pf", "c_pd",
                 "c_subseq", "c_terminal", "ae_cost")

# One evaluation with modified spec/HRs
evaluate <- function(pieces, spec, hr = NULL) {
  cv <- assemble_curves(pieces, spec, hr)
  r <- run_engine(cv, spec)
  list(curves = cv, res = r)
}

# ---- deterministic one-way sensitivity analysis -------------------------------
dsa_params <- function(spec, cost_pct = 0.2, wane_delta = 2) {
  p <- list()
  add <- function(name, low, high, setter) {
    p[[length(p) + 1]] <<- list(name = name, low = low, high = high, set = setter)
  }
  for (ep in c("OS", "PFS")) {
    h <- spec$endpoints[[ep]]$hr
    if (spec$endpoints[[ep]]$approach == "hr" && is.finite(h$se)) {
      add(paste("HR", ep), h$lower, h$upper,
          local({ e <- ep; \(s, v) { s$endpoints[[e]]$hr$hr <- v; s } }))
    }
  }
  u <- spec$econ$utility
  add("Utility, progression-free", max(u$pf - 1.96 * u$pf_se, 0), min(u$pf + 1.96 * u$pf_se, 1),
      \(s, v) { s$econ$utility$pf <- v; s })
  add("Utility, progressed", max(u$pd - 1.96 * u$pd_se, 0), min(u$pd + 1.96 * u$pd_se, 1),
      \(s, v) { s$econ$utility$pd <- v; s })
  labels <- c(upfront = "upfront cost", drug_cost = "drug cost", admin_cost = "admin cost",
              c_pf = "PF care cost", c_pd = "PD care cost", c_subseq = "subsequent therapy",
              c_terminal = "terminal care", ae_cost = "AE cost")
  for (a in c("int", "comp")) for (f in cost_fields) {
    base <- spec$econ[[a]][[f]]
    if (is.finite(base) && base > 0) {
      add(paste0(spec$arms[[a]], ": ", labels[[f]]), base * (1 - cost_pct), base * (1 + cost_pct),
          local({ aa <- a; ff <- f; \(s, v) { s$econ[[aa]][[ff]] <- v; s } }))
    }
  }
  # Waning timing (V3): structural, but reviewers test its timing (Doc 12, TA798, TA1099)
  wn <- spec$waning
  if (wn$type != "none" && is.finite(wn$start) && wane_delta > 0) {
    add("Waning start (years, end moves with it)", max(0, wn$start - wane_delta), wn$start + wane_delta,
        \(s, v) { d <- s$waning$end - s$waning$start; s$waning$start <- v; s$waning$end <- v + d; s })
    if (wn$type == "gradual") {
      dur <- wn$end - wn$start
      add("Waning duration (years)", max(0.25, dur - wane_delta), dur + wane_delta,
          \(s, v) { s$waning$end <- s$waning$start + v; s })
    }
  }
  add("Discount rate, costs", 0, 2 * spec$settings$dr_cost,
      \(s, v) { s$settings$dr_cost <- v; s })
  add("Discount rate, outcomes", 0, 2 * spec$settings$dr_out,
      \(s, v) { s$settings$dr_out <- v; s })
  p
}

run_dsa <- function(pieces, spec, cost_pct = 0.2, wane_delta = 2) {
  base <- evaluate(pieces, spec)$res
  rows <- lapply(dsa_params(spec, cost_pct, wane_delta), \(p) {
    one <- \(v) {
      s2 <- p$set(spec, v)
      hr <- list(OS = s2$endpoints$OS$hr$hr, PFS = s2$endpoints$PFS$hr$hr)
      evaluate(pieces, s2, hr)$res
    }
    lo <- one(p$low)
    hi <- one(p$high)
    data.frame(parameter = p$name, low_value = p$low, high_value = p$high,
               icer_low = lo$icer, icer_high = hi$icer,
               nmb_low = lo$nmb, nmb_high = hi$nmb)
  })
  out <- do.call(rbind, rows)
  attr(out, "base") <- base
  out
}

# ---- probabilistic sensitivity analysis -----------------------------------------
# Sampled: log-HRs (normal, optional OS/PFS correlation), utilities (beta),
# costs (gamma) and, new in V2, the fitted survival curves (surv_method):
#   "none" curves fixed (V1 behaviour) | "mvn" parameter covariance | "boot" bootstrap refit
# Parameter draws come first, so for surv_method = "none" V2 reproduces V1 exactly.
run_psa <- function(pieces, spec, n = 1000, seed = 1, cost_cv = 0.2, hr_corr = 0,
                    surv_method = c("none", "mvn", "boot"), fits = NULL, ipd = NULL,
                    cores = 1, code_dir = NULL, band_step = 0.25, batch = 25,
                    wane_start_range = NULL, progress = NULL) {
  surv_method <- match.arg(surv_method)
  set.seed(seed)
  h_os <- spec$endpoints$OS$hr
  h_pfs <- spec$endpoints$PFS$hr
  use_os <- spec$endpoints$OS$approach == "hr" && is.finite(h_os$se)
  use_pfs <- spec$endpoints$PFS$approach == "hr" && is.finite(h_pfs$se)
  z1 <- stats::rnorm(n)
  z2 <- hr_corr * z1 + sqrt(1 - hr_corr^2) * stats::rnorm(n)
  u <- spec$econ$utility
  draws <- list(
    hr_os = if (use_os) exp(log(h_os$hr) + h_os$se * z1) else rep(h_os$hr %||% 1, n),
    hr_pfs = if (use_pfs) exp(log(h_pfs$hr) + h_pfs$se * z2) else rep(h_pfs$hr %||% 1, n),
    u_pf = rbeta_ms(n, u$pf, u$pf_se),
    u_pd = rbeta_ms(n, u$pd, u$pd_se),
    costs = lapply(c(comp = "comp", int = "int"), \(a) {
      sapply(cost_fields, \(f) rgamma_mcv(n, spec$econ[[a]][[f]], cost_cv))
    })
  )
  draws$run_seed <- sample.int(.Machine$integer.max, n)
  sampler <- make_surv_sampler(spec, fits, pieces, surv_method, n, ipd)
  # V3: optional uniform draw of the waning start (duration kept). Drawn last, so
  # every other draw is identical to V2 for the same seed.
  if (!is.null(wane_start_range) && spec$waning$type != "none") {
    draws$wane_start <- stats::runif(n, min(wane_start_range), max(wane_start_range))
  }

  band_t <- seq(0, spec$settings$horizon, by = band_step)
  band_idx <- pmax(findInterval(band_t, pieces$t), 1)
  run_fun <- make_run_fun(spec, sampler, draws, band_idx)

  cores <- max(1, as.integer(cores))
  if (cores > 1) {
    if (is.null(code_dir)) stop("Parallel PSA needs code_dir (the app's R/ folder).")
    cl <- start_psa_cluster(cores, code_dir, run_fun)
    on.exit(parallel::stopCluster(cl), add = TRUE)
    worker <- eval(quote(function(i) .psm_run(i)), globalenv())
  }
  out <- vector("list", n)
  done <- 0
  for (b in split(seq_len(n), ceiling(seq_len(n) / batch))) {
    out[b] <- if (cores > 1) parallel::parLapply(cl, b, worker) else lapply(b, run_fun)
    done <- done + length(b)
    if (!is.null(progress)) progress(done / n)
  }

  runs <- do.call(rbind, lapply(out, `[[`, "row"))
  runs$d_cost <- runs$cost_int - runs$cost_comp
  runs$d_qaly <- runs$qaly_int - runs$qaly_comp
  bands <- do.call(rbind, lapply(out, `[[`, "band"))
  k <- length(band_t)
  attr(runs, "bands") <- list(
    t = band_t,
    curves = list(comp_OS = bands[, 1:k, drop = FALSE], comp_PFS = bands[, k + 1:k, drop = FALSE],
                  int_OS = bands[, 2 * k + 1:k, drop = FALSE], int_PFS = bands[, 3 * k + 1:k, drop = FALSE])
  )
  attr(runs, "settings") <- list(n = n, seed = seed, cost_cv = cost_cv, hr_corr = hr_corr,
                                 surv_method = surv_method, cores = cores,
                                 wane_start_range = if (!is.null(draws$wane_start)) range(wane_start_range) else NULL,
                                 paired = !is.null(ipd) && ipd_is_paired(ipd[ipd$arm %in% unlist(spec$arms), ]),
                                 sampled_curves = vapply(surv_tasks(spec), \(x) paste(x$ep, paste(spec$arms[x$arms], collapse = "+")), ""))
  runs
}

# Median and 95% interval of the sampled curves (valid runs only)
psa_bands <- function(runs, arms) {
  b <- attr(runs, "bands")
  if (is.null(b)) return(NULL)
  ok <- runs$valid
  lab <- c(comp_OS = "OS", comp_PFS = "PFS", int_OS = "OS", int_PFS = "PFS")
  arm <- c(comp_OS = arms$comp, comp_PFS = arms$comp, int_OS = arms$int, int_PFS = arms$int)
  do.call(rbind, lapply(names(b$curves), \(nm) {
    m <- b$curves[[nm]][ok, , drop = FALSE]
    q <- apply(m, 2, stats::quantile, probs = c(0.025, 0.5, 0.975), na.rm = TRUE)
    data.frame(time = b$t, arm = arm[[nm]], endpoint = lab[[nm]],
               lower = q[1, ], median = q[2, ], upper = q[3, ])
  }))
}

# Deterministic result vs the mean of valid PSA runs
psa_vs_deterministic <- function(runs, res) {
  ok <- runs[runs$valid, ]
  det <- c(cost_comp = res$comp$summary[["cost_total"]], cost_int = res$int$summary[["cost_total"]],
           qaly_comp = res$comp$summary[["qaly"]], qaly_int = res$int$summary[["qaly"]],
           ly_comp = res$comp$summary[["ly_disc"]], ly_int = res$int$summary[["ly_disc"]],
           d_cost = res$inc_cost, d_qaly = res$inc_qaly)
  prob <- vapply(names(det), \(k) mean(ok[[k]]), numeric(1))
  data.frame(measure = names(det), deterministic = det, probabilistic = prob,
             rel_diff = ifelse(abs(det) > 1e-12, prob / det - 1, NA), row.names = NULL)
}

# Running mean of increments over valid runs, in run order (convergence check)
psa_convergence <- function(runs, wtp) {
  ok <- runs[runs$valid, ]
  k <- seq_len(nrow(ok))
  mc <- cumsum(ok$d_cost) / k
  mq <- cumsum(ok$d_qaly) / k
  data.frame(runs = k, icer = ifelse(abs(mq) > 1e-12, mc / mq, NA), nmb = wtp * mq - mc)
}

# Summary: ICER as the ratio of mean increments (Doc 13 Case 17)
psa_summary <- function(runs, wtp) {
  ok <- runs[runs$valid, ]
  m_cost <- if (nrow(ok)) mean(ok$d_cost) else NA_real_
  m_qaly <- if (nrow(ok)) mean(ok$d_qaly) else NA_real_
  list(
    n = nrow(runs), n_valid = nrow(ok), n_invalid = sum(!runs$valid),
    n_constraint = sum(runs$constraint_active),
    mean_d_cost = m_cost, mean_d_qaly = m_qaly,
    icer = if (is.finite(m_qaly) && abs(m_qaly) > 1e-12) m_cost / m_qaly else NA,
    p_ce = if (nrow(ok)) mean(wtp * ok$d_qaly - ok$d_cost > 0) else NA_real_,
    reasons = if (any(!runs$valid)) table(runs$reason[!runs$valid]) else NULL,
    settings = attr(runs, "settings")
  )
}

ceac <- function(runs, wtp_max, steps = 60) {
  ok <- runs[runs$valid, ]
  w <- seq(0, wtp_max, length.out = steps)
  data.frame(wtp = w, p_ce = vapply(w, \(x) mean(x * ok$d_qaly - ok$d_cost > 0), numeric(1)))
}
