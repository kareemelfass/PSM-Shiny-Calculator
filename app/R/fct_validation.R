# fct_validation.R : automated checks ---------------------------------------
# Domains follow Doc 13 (structural, input, computational).
# Status: pass / warn / fail / info. Nothing is fixed silently: constraints the
# model applies are reported here with when and where they bind.

check_row <- function(domain, check, status, detail, source = "") {
  data.frame(domain = domain, check = check, status = status,
             detail = detail, source = source)
}

curve_checks <- function(curves, arms) {
  rows <- list()
  bad_range <- character()
  bad_mono <- character()
  for (a in c("comp", "int")) for (ep in c("OS", "PFS", "ToT")) {
    s <- curves[[a]][[ep]]
    lab <- paste(arms[[a]], ep)
    if (any(!is.finite(s)) || any(s < -1e-12 | s > 1 + 1e-12)) bad_range <- c(bad_range, lab)
    if (any(diff(s) > 1e-10)) bad_mono <- c(bad_mono, lab)
  }
  rows[[1]] <- check_row("Computational", "Curves within 0 to 1",
    if (length(bad_range)) "fail" else "pass",
    if (length(bad_range)) paste("Out of range:", paste(bad_range, collapse = ", ")) else
      "All OS, PFS and ToT curves stay within 0 to 1 over the full horizon.",
    "Doc 3; Doc 13 Case 16")
  rows[[2]] <- check_row("Computational", "Survival never increases (full horizon)",
    if (length(bad_mono)) "fail" else "pass",
    if (length(bad_mono)) paste("Increasing at some cycle:", paste(bad_mono, collapse = ", ")) else
      "Checked every cycle to the end of the horizon, not just the observed period.",
    "Doc 13 Case 14")
  do.call(rbind, rows)
}

occupancy_ok <- function(curves, tol = 1e-9) {
  all(vapply(c("comp", "int"), \(a) {
    cv <- curves[[a]]
    pf <- cv$PFS; pd <- cv$OS - cv$PFS; dead <- 1 - cv$OS
    all(pf >= -tol, pd >= -tol, dead >= -tol, abs(pf + pd + dead - 1) < tol)
  }, logical(1)))
}

validate_model <- function(spec, pieces, curves, res, fits = NULL) {
  arms <- spec$arms
  st <- spec$settings
  out <- list(curve_checks(curves, arms))
  add <- function(...) out[[length(out) + 1]] <<- check_row(...)

  # Occupancy
  add("Computational", "State occupancy valid and sums to 1",
      if (occupancy_ok(curves)) "pass" else "fail",
      "PF = PFS, PD = OS - PFS, Dead = 1 - OS; each between 0 and 1, sum = 1, every cycle, both arms.",
      "Doc 8; Doc 13 Case 16")

  # HR = 1 must reproduce the reference curve (live engine test)
  test_spec <- spec
  for (ep in c("OS", "PFS")) test_spec$endpoints[[ep]]$approach <- "hr"
  test_spec$waning$type <- "none"
  test_spec$settings$floor_on <- FALSE
  tp <- pieces
  tp$gp <- NULL
  tc <- assemble_curves(tp, test_spec, hr = list(OS = 1, PFS = 1))
  diff_max <- max(abs(tc$int$OS - tc$comp$OS), abs(tc$int$PFS - tc$comp$PFS))
  add("Computational", "HR = 1 reproduces the reference curve",
      if (diff_max < 1e-12) "pass" else "fail",
      paste0("Maximum difference: ", signif(diff_max, 3), "."), "Doc 3")

  # Totals add up; QALYs and discounting sensible
  for (a in c("comp", "int")) {
    sm <- res[[a]]$summary
    parts <- sum(sm[grep("^cost_(?!total)", names(sm), perl = TRUE)])
    ok <- abs(parts - sm[["cost_total"]]) < 1e-6 &&
      sm[["qaly"]] <= sm[["ly_disc"]] + 1e-9 &&
      sm[["ly_disc"]] <= sm[["ly_undisc"]] + 1e-9
    add("Computational", paste("Totals consistent:", arms[[a]]),
        if (ok) "pass" else "fail",
        "Cost components sum to the total; QALYs do not exceed life-years; discounted values do not exceed undiscounted.",
        "Doc 13")
  }

  # PFS above OS before the constraint
  for (a in c("comp", "int")) {
    f <- curves$flags[[a]]
    n_x <- f$pfs_gt_os_cycles %||% 0
    add("Structural", paste("PFS not above OS:", arms[[a]]),
        if (n_x == 0) "pass" else if (isTRUE(st$pfs_le_os)) "warn" else "fail",
        if (n_x == 0) "Fitted PFS stays at or below OS in every cycle." else
          paste0("Fitted PFS exceeds OS from year ", round(f$pfs_gt_os_from, 2), " (", n_x,
                 " cycles). ", if (isTRUE(st$pfs_le_os))
                   "The explicit PFS <= OS constraint is applied there. Consider jointly plausible OS/PFS models." else
                   "No constraint applied: occupancy will be invalid."),
        "Doc 9 C3; Doc 13 Case 16")
    n_t <- f$tot_gt_os_cycles %||% 0
    add("Structural", paste("ToT not above OS:", arms[[a]]),
        if (n_t == 0) "pass" else "warn",
        if (n_t == 0) "Time on treatment never exceeds overall survival." else
          paste0("Fitted ToT exceeds OS from year ", round(f$tot_gt_os_from, 2),
                 "; capped at OS (no one is treated after death)."),
        "Doc 11")
    if (!isTRUE(spec$tot[[a]]$beyond_prog) && (f$tot_gt_pfs_cycles %||% 0) > 0) {
      add("Structural", paste("ToT capped at PFS:", arms[[a]]), "info",
          paste0("Treatment beyond progression is off, so ToT was capped at PFS in ",
                 f$tot_gt_pfs_cycles, " cycles."), "Doc 11")
    }
  }

  # Crossover coherence (Doc 10)
  cx <- spec$crossover
  if (!is.null(cx)) {
    incoherent <- (cx$comp_os_basis == "itt" && !isTRUE(cx$subseq_includes_int)) ||
      (cx$comp_os_basis == "adjusted" && isTRUE(cx$subseq_includes_int))
    add("Structural", "Crossover: survival and costs describe the same pathway",
        if (incoherent) "warn" else "pass",
        if (!incoherent) "Comparator survival basis and subsequent-therapy costs are aligned." else
          if (cx$comp_os_basis == "itt") "ITT comparator OS includes benefit from switching to the intervention, but its cost is not included in comparator subsequent therapy." else
            "Comparator OS is adjusted to a no-switching pathway, yet comparator costs still include the intervention after switching.",
        "Doc 10")
  }

  # Treatment-effect duration and waning (Docs 11, 12; V3)
  wn <- spec$waning
  cap <- spec$tot$int$cap %||% NA
  dg <- waning_diagnostics(curves, spec)
  fu <- st$fu_years %||% NA
  if (wn$type == "none") {
    add("Input", "Sustained treatment effect",
        if (is.finite(cap) && cap > 0) "warn" else "info",
        if (is.finite(cap) && cap > 0) paste0("Intervention treatment stops by year ", cap,
          " but its relative effect is assumed to last the whole horizon. A stopping rule does not justify",
          " lifetime efficacy (NICE TA650): provide evidence of durable benefit or test waning scenarios.") else
          "The relative effect lasts the whole horizon. This needs evidence of durable benefit (long-term data, cure, mechanism).",
        "Doc 12; Doc 11")
  } else {
    lab <- waning_structures()$label[waning_structures()$type == wn$type]
    timing <- if (wn$type == "gradual") sprintf("from year %s to year %s", round(dg$start, 2), round(dg$end, 2)) else
      sprintf("at year %s", round(dg$start, 2))
    add("Input", "Waning assumption", if (wn$type == "stop_drop") "warn" else "info",
        paste0(lab, " ", timing, if (isTRUE(dg$anchored_to_cap)) paste0(" (", wn$start, " years after treatment stops at year ", cap, ")") else
          " (years since treatment start)", ", applied to ",
          paste(c("OS", "PFS")[c(isTRUE(wn$apply_os), isTRUE(wn$apply_pfs))], collapse = " and "),
          if ((wn$target_hr %||% 1) != 1) paste0("; HR after waning = ", wn$target_hr) else "",
          if ((wn$frailty %||% 0) > 0) paste0("; individual-level waning with frailty variance ", wn$frailty) else "", ". ",
          waning_structures()$role[waning_structures()$type == wn$type], "."),
        "Doc 12")
    if (isTRUE(dg$anchor_missing_cap)) {
      add("Input", "Waning anchored to treatment stop", "warn",
          "Timing is set relative to treatment stop, but no maximum treatment duration is set on the Time on treatment tab. Timing was measured from treatment start instead.",
          "Doc 12")
    }
    if (is.finite(fu) && is.finite(dg$start) && dg$start < fu) {
      add("Input", "Waning starts inside the observed period", "warn",
          sprintf("Waning starts at year %s, before trial follow-up ends (year %s). The model then overrides observed data; waning is an assumption about the unobserved period.",
                  round(dg$start, 2), round(fu, 2)),
          "Doc 12")
    }
    for (ep in c("OS", "PFS")) {
      d <- dg[[ep]]
      if (is.null(d)) next
      if (!is.null(d$drop_pts)) {
        add("Structural", paste("Stop and drop: instant loss,", ep), "warn",
            sprintf("At year %s, %s points of %s (%s of intervention patients still alive) are lost instantly. This is a survival jump, not a finite hazard: extreme stress test only.",
                    round(dg$start, 2), fmt_num(100 * d$drop_pts, 1), ep, fmt_pct(d$drop_share)),
            "Doc 12")
      }
      if (!is.null(d$hr_after) && (wn$frailty %||% 0) > 0) {
        add("Input", paste("Marginal HR after individual-level waning,", ep), "info",
            sprintf("Population-level HR just after waning = %.2f (target %s). Survivors in the intervention arm are frailer, so equal individual hazards give a higher population hazard. Constraining the marginal HR to 1 would overestimate benefit (Jennings 2024).",
                    d$hr_after, wn$target_hr %||% 1),
            "Doc 12 ref 2")
      }
      if (!is.null(d$share_above_target) && d$share_above_target > 0) {
        material <- d$share_above_weighted > 0.05
        add("Input", paste("Waning lowers the intervention's hazard,", ep),
            if (material) "warn" else "info",
            sprintf("From year %s the fitted intervention hazard is above the target, so waning lowers it (favours the intervention). %s of intervention patients remain at that point; weighted by patients at risk this affects %s of the waning period%s",
                    round(d$first_above, 1), fmt_pct(d$at_risk_first_above), fmt_pct(d$share_above_weighted),
                    if (material) ". Check the independent fits on the implied-HR plot." else ", so the impact is small."),
            "Doc 12; Doc 3")
      }
    }
  }
  if (!nzchar(trimws(wn$rationale %||% ""))) {
    add("Input", "Justification for the effect-duration assumption", "warn",
        "None recorded. HTA committees judge the justification, not the method: record the clinical evidence and the scenarios tested (Treatment effects tab).",
        "Doc 12")
  } else {
    add("Input", "Justification for the effect-duration assumption", "pass",
        paste0("Recorded: ", substr(wn$rationale, 1, 300)), "Doc 12")
  }

  # HR direction statement
  for (ep in c("OS", "PFS")) {
    es <- spec$endpoints[[ep]]
    if (es$approach == "hr") {
      h <- es$hr
      add("Input", paste("HR applied to", ep), "info",
          paste0("HR intervention vs comparator = ", round(h$hr, 3), " (",
                 if (h$hr < 1) paste0(round(100 * (1 - h$hr), 1), "% lower hazard") else
                   paste0(round(100 * (h$hr - 1), 1), "% higher hazard"),
                 "), applied as S_int(t) = S_comp(t)^HR to the ", ep, " curve only."),
          "Doc 3; Doc 9 B4, C2")
    }
  }

  # General-population mortality floor
  if (isTRUE(st$floor_on) && !is.null(curves$gp)) {
    det <- vapply(c("comp", "int"), \(a) {
      f <- curves$flags[[a]]$floor_from
      paste0(arms[[a]], ": ", if (is.na(f)) "never binds" else paste("binds from year", round(f, 1)))
    }, character(1))
    add("Input", "General-population mortality floor", "info", paste(det, collapse = "; "), "Doc 7")
  } else {
    add("Input", "General-population mortality floor", "warn",
        "Off: modelled mortality may fall below that of the general population.", "Doc 7")
  }

  # Hybrid cut-point risk set (Doc 6)
  if (!is.null(fits)) for (ep in c("OS", "PFS")) {
    es <- spec$endpoints[[ep]]
    if (isTRUE(es$hybrid)) for (a in names(fits[[ep]]$km)) {
      if (es$approach == "hr" && a == "int") next
      km <- fits[[ep]]$km[[a]]
      n0 <- km$n
      nr <- km_at_risk(km, es$cut)
      add("Input", paste("Hybrid cut-point risk set:", ep, arms[[a]]),
          if (nr < 20 || nr < 0.1 * n0) "warn" else "pass",
          paste0(nr, " of ", n0, " patients at risk at year ", es$cut,
                 ". Justify the cut-point and test alternatives."),
          "Doc 6")
    }
  }

  # Extrapolation dependence (Doc 9 C4a)
  sh <- res$share_ly_after_fu
  if (is.finite(sh)) {
    add("Input", "Incremental life-years gained after trial follow-up",
        if (sh > 0.5) "warn" else "info",
        paste0(fmt_pct(sh), " of incremental (undiscounted) life-years accrue after year ",
               round(st$fu_years, 2), ", outside the observed data."),
        "Doc 9 C4a; Doc 6")
  }

  # PH evidence where a single HR or joint PH-type model is used
  if (!is.null(fits)) for (ep in c("OS", "PFS")) {
    es <- spec$endpoints[[ep]]
    if (es$approach %in% c("hr", "joint") && all(fits[[ep]]$have)) {
      cx <- tryCatch(cox_summary(fits[[ep]]), error = \(e) NULL)
      if (!is.null(cx) && is.finite(cx$p_ph)) {
        add("Input", paste("Proportional hazards,", ep),
            if (cx$p_ph < 0.05) "warn" else "info",
            paste0("Schoenfeld test p = ", signif(cx$p_ph, 3),
                   if (cx$p_ph < 0.05) ": evidence against PH." else
                     ". A non-significant test does not prove PH; inspect the plots."),
            "Doc 3")
      }
    }
  }

  do.call(rbind, out)
}

# Landmark survival vs clinical benchmark ranges
landmark_table <- function(curves, arms, years, bench = NULL) {
  t <- curves$t
  at <- \(s, y) stats::approx(t, s, xout = y, rule = 2)$y
  rows <- lapply(c("comp", "int"), \(a) lapply(c("OS", "PFS"), \(ep) {
    v <- at(curves[[a]][[ep]], years)
    data.frame(arm = arms[[a]], endpoint = ep, t(setNames(v, paste0("yr_", years))),
               check.names = FALSE)
  }))
  tb <- do.call(rbind, unlist(rows, recursive = FALSE))
  tb
}

parse_bench <- function(txt, n) {
  if (is.null(txt) || !nzchar(trimws(txt))) return(NULL)
  parts <- trimws(strsplit(txt, ",")[[1]])
  if (length(parts) != n) return(NULL)
  rng <- lapply(parts, \(p) {
    v <- suppressWarnings(as.numeric(strsplit(p, "-")[[1]]))
    if (length(v) != 2 || any(!is.finite(v))) c(NA, NA) else v / 100
  })
  do.call(rbind, rng)
}
