# fct_psa_surv.R : survival-curve uncertainty for the PSA (new in V2) ---------
# A "sampler" returns, for PSA run i, a full set of model pieces (OS, PFS, ToT
# curves on the model grid) drawn from their uncertainty. Two methods:
#   mvn  : draw model parameters from N(estimate, covariance). Fast.
#   boot : resample patients with replacement, refit the chosen models. Slower;
#          with patient ids, OS, PFS and ToT are resampled together, which keeps
#          their correlation (Doc 13 Case 16).
# Only curves that come from a fitted model are sampled. Hazard ratios entered
# by the user are sampled separately in run_psa().

# Which fitted curves does this specification use?
surv_tasks <- function(spec) {
  tasks <- list()
  add <- \(...) tasks[[length(tasks) + 1]] <<- list(...)
  for (ep in c("OS", "PFS")) {
    es <- spec$endpoints[[ep]]
    if (es$approach == "joint") {
      add(ep = ep, arms = c("comp", "int"), joint = TRUE, model = es$comp_model,
          hybrid = isTRUE(es$hybrid), cut = es$cut)
    } else {
      add(ep = ep, arms = "comp", joint = FALSE, model = es$comp_model,
          hybrid = isTRUE(es$hybrid), cut = es$cut)
      if (es$approach == "indep") {
        add(ep = ep, arms = "int", joint = FALSE, model = es$int_model,
            hybrid = isTRUE(es$hybrid), cut = es$cut)
      }
    }
  }
  for (a in c("comp", "int")) {
    ts <- spec$tot[[a]]
    if (identical(ts$method, "fit")) {
      add(ep = "ToT", arms = a, joint = FALSE, model = ts$model, hybrid = FALSE, cut = NA)
    }
  }
  tasks
}

put_piece <- function(p, ep, a, s) {
  p[[ep]][[a]] <- s
  p
}

# ---- bootstrap resampling -----------------------------------------------------

# Resample within arm. Rows with a patient id are resampled as whole patients;
# rows without one are resampled within arm x endpoint.
boot_ipd <- function(ipd) {
  if (!"id" %in% names(ipd)) ipd$id <- NA_character_
  has_id <- !is.na(ipd$id)
  idx <- integer()
  paired <- which(has_id)
  if (length(paired)) {
    by_arm <- split(paired, ipd$arm[paired])
    for (rows in by_arm) {
      by_id <- split(rows, ipd$id[rows])
      pick <- sample(names(by_id), length(by_id), replace = TRUE)
      idx <- c(idx, unlist(by_id[pick], use.names = FALSE))
    }
  }
  loose <- which(!has_id)
  if (length(loose)) {
    strata <- split(loose, paste(ipd$arm[loose], ipd$endpoint[loose]))
    for (rows in strata) idx <- c(idx, rows[sample.int(length(rows), length(rows), replace = TRUE)])
  }
  ipd[idx, ]
}

# ---- the sampler ---------------------------------------------------------------

make_surv_sampler <- function(spec, fits, pieces, method = c("none", "mvn", "boot"), n = 1,
                              ipd = NULL) {
  # Force every argument: the returned closure is shipped to parallel workers,
  # and an unevaluated argument (a promise) cannot be resolved there.
  force(spec); force(fits); force(pieces); force(n); force(ipd)
  method <- match.arg(method)
  tasks <- surv_tasks(spec)
  if (method == "none" || !length(tasks)) {
    return(function(i, seed = NULL) list(pieces = pieces, fail = NULL))
  }
  t <- pieces$t
  arms <- spec$arms
  fit_for <- function(task, a) {
    fe <- fits[[task$ep]]
    f <- if (task$joint) fe$joint[[task$model]] else fe$indep[[a]][[task$model]]
    if (is.null(f) || !isTRUE(f$ok)) stop(task$ep, " model ", task$model, " is not available")
    f
  }

  if (method == "mvn") {
    # Draw all parameters now: runs are then identical in serial and parallel
    draws <- list()
    for (k in seq_along(tasks)) {
      task <- tasks[[k]]
      if (task$joint) {
        draws[[k]] <- list(all = fit_for(task, "comp")$draw_mvn(n))
      } else {
        a <- task$arms
        draws[[k]] <- setNames(list(fit_for(task, a)$draw_mvn(n)), a)
      }
    }
    return(function(i, seed = NULL) {
      p <- pieces
      for (k in seq_along(tasks)) {
        task <- tasks[[k]]
        for (a in task$arms) {
          dfun <- if (task$joint) draws[[k]]$all else draws[[k]][[a]]
          lab <- arms[[a]]
          s <- hybrid_curve(\(x) dfun(i, x, lab), t, fits[[task$ep]]$km[[a]], task$hybrid, task$cut)
          p <- put_piece(p, task$ep, a, s)
        }
      }
      list(pieces = p, fail = NULL)
    })
  }

  # Bootstrap
  if (is.null(ipd)) stop("The bootstrap needs the patient-level data.")
  ipd <- ipd[ipd$arm %in% unlist(arms), ]
  levels_arm <- c(arms$comp, arms$int)
  function(i, seed = NULL) {
    if (!is.null(seed)) set.seed(seed)
    b <- boot_ipd(ipd)
    p <- pieces
    for (task in tasks) {
      engine <- fits[[task$ep]]$engine
      d_all <- b[b$endpoint == task$ep, c("arm", "time", "event")]
      d_all$arm <- factor(d_all$arm, levels = levels_arm)
      joint_fit <- if (task$joint) fit_one(d_all, task$model, TRUE, engine) else NULL
      for (a in task$arms) {
        d <- d_all[d_all$arm == arms[[a]], ]
        d$arm <- droplevels(d$arm)
        f <- joint_fit %||% fit_one(d, task$model, FALSE, engine)
        if (!isTRUE(f$ok)) {
          return(list(pieces = NULL, fail = paste("Refit failed:", task$ep, arms[[a]], model_label(task$model))))
        }
        km <- if (task$hybrid) survival::survfit(survival::Surv(time, event) ~ 1, data = d) else NULL
        lab <- arms[[a]]
        s <- hybrid_curve(\(x) f$S(x, lab), t, km, task$hybrid, task$cut)
        p <- put_piece(p, task$ep, a, s)
      }
    }
    list(pieces = p, fail = NULL)
  }
}

# ---- one PSA run ----------------------------------------------------------------

# Built in its own function so the closure carries only what a run needs
# (it is shipped to parallel workers; it must not capture a Shiny session).
make_run_fun <- function(spec, sampler, draws, band_idx) {
  force(spec); force(sampler); force(draws); force(band_idx)
  function(i) {
    sp <- tryCatch(sampler(i, draws$run_seed[i]),
                   error = \(e) list(pieces = NULL, fail = paste("Sampling error:", conditionMessage(e))))
    row <- data.frame(run = i, hr_os = draws$hr_os[i], hr_pfs = draws$hr_pfs[i],
                      wane_start = if (is.null(draws$wane_start)) NA_real_ else draws$wane_start[i],
                      u_pf = draws$u_pf[i], u_pd = draws$u_pd[i],
                      cost_comp = NA_real_, cost_int = NA_real_, qaly_comp = NA_real_,
                      qaly_int = NA_real_, ly_comp = NA_real_, ly_int = NA_real_,
                      valid = FALSE, constraint_active = FALSE, reason = "")
    if (!is.null(sp$fail)) {
      row$reason <- sp$fail
      return(list(row = row, band = rep(NA_real_, 4 * length(band_idx))))
    }
    s2 <- spec
    s2$econ$utility$pf <- draws$u_pf[i]
    s2$econ$utility$pd <- draws$u_pd[i]
    if (!is.null(draws$wane_start)) {
      d <- s2$waning$end - s2$waning$start
      s2$waning$start <- draws$wane_start[i]
      s2$waning$end <- draws$wane_start[i] + d
    }
    for (a in c("comp", "int")) for (f in cost_fields) s2$econ[[a]][[f]] <- draws$costs[[a]][i, f]
    ev <- evaluate(sp$pieces, s2, list(OS = draws$hr_os[i], PFS = draws$hr_pfs[i]))
    cv <- ev$curves
    fl <- cv$flags
    rising <- vapply(c("comp", "int"), \(a) {
      any(diff(cv[[a]]$OS) > 1e-10) || any(diff(cv[[a]]$PFS) > 1e-10)
    }, logical(1))
    occ <- occupancy_ok(cv)
    row$cost_comp <- ev$res$comp$summary[["cost_total"]]
    row$cost_int <- ev$res$int$summary[["cost_total"]]
    row$qaly_comp <- ev$res$comp$summary[["qaly"]]
    row$qaly_int <- ev$res$int$summary[["qaly"]]
    row$ly_comp <- ev$res$comp$summary[["ly_disc"]]
    row$ly_int <- ev$res$int$summary[["ly_disc"]]
    row$valid <- occ && !any(rising) && is.finite(row$cost_int + row$qaly_int)
    row$constraint_active <- (fl$comp$pfs_gt_os_cycles %||% 0) + (fl$int$pfs_gt_os_cycles %||% 0) > 0
    row$reason <- if (row$valid) "" else if (any(rising)) {
      paste("Survival increases:", paste(spec$arms[names(rising)[rising]], collapse = ", "))
    } else if (!occ) "State occupancy invalid" else "Non-finite result"
    band <- c(cv$comp$OS[band_idx], cv$comp$PFS[band_idx], cv$int$OS[band_idx], cv$int$PFS[band_idx])
    list(row = row, band = band)
  }
}

# ---- parallel workers (base R 'parallel'; works on Windows) ---------------------

start_psa_cluster <- function(cores, code_dir, run_fun) {
  cl <- parallel::makePSOCKcluster(cores)
  ok <- FALSE
  on.exit(if (!ok) parallel::stopCluster(cl))
  load_code <- function(dir) {
    files <- sort(list.files(dir, pattern = "^(fct_|utils).*\\.[Rr]$", full.names = TRUE))
    for (f in files) sys.source(f, envir = globalenv())
    requireNamespace("survival", quietly = TRUE)
    NULL
  }
  set_run <- function(f) {
    assign(".psm_run", f, envir = globalenv())
    NULL
  }
  environment(load_code) <- globalenv()   # ship the bare function, not this frame
  environment(set_run) <- globalenv()
  parallel::clusterCall(cl, load_code, code_dir)
  parallel::clusterCall(cl, set_run, run_fun)
  ok <- TRUE
  cl
}
