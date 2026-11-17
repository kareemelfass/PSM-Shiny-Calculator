# fct_waning.R : waning diagnostics and the waning explorer (V3, Doc 12) -----------

waning_structures <- function() {
  data.frame(
    type = c("none", "equal", "gradual", "stop_drop"),
    label = c("Sustained effect", "Equal hazards from start", "Gradual equalisation",
              "Stop and drop"),
    role = c("Base case only with evidence of durable benefit; otherwise a scenario",
             "Pragmatic base-case candidate; test timing in scenarios",
             "Reasonable base-case candidate; test timing and duration",
             "Extreme stress test, not a base case")
  )
}

# Diagnostics for the current waning setting (used by plots and validation)
waning_diagnostics <- function(curves, spec) {
  wn <- curves$waning
  t <- curves$t
  out <- list(type = wn$type, start = wn$start, end = wn$end,
              anchored_to_cap = wn$anchored_to_cap, anchor_missing_cap = wn$anchor_missing_cap)
  if (wn$type == "none" || !is.finite(wn$start)) return(out)
  tm <- cycle_mid(t)
  for (ep in c("OS", "PFS")) {
    if (!isTRUE(wn[[paste0("apply_", tolower(ep))]])) next
    s_c <- curves$comp[[ep]]
    s_pre <- curves$pre_waning[[ep]]
    s_int <- curves$int[[ep]]
    d <- list()
    # Stop and drop: instantaneous loss at the start time (not a finite hazard)
    if (wn$type == "stop_drop") {
      i <- which(t >= wn$start)[1]
      if (!is.na(i) && i > 1) {
        d$drop_pts <- s_pre[i - 1] - s_c[i - 1]
        d$drop_share <- d$drop_pts / max(s_pre[i - 1], 1e-12)
      }
    }
    # Marginal HR implied just after waning is complete
    after <- which(tm >= max(wn$end, wn$start))[1:4]
    after <- after[!is.na(after)]
    if (length(after)) {
      hi <- cycle_hazard(s_int)[after]
      hc <- cycle_hazard(s_c)[after]
      d$hr_after <- mean(hi / pmax(hc, 1e-300))
    }
    # Direction: before waning, is the intervention hazard above the target?
    w <- waning_weight(tm, wn$start, if (wn$type == "gradual") wn$end else wn$start)
    zone <- which(w > 0)
    if (length(zone)) {
      raw <- cycle_hazard(s_pre)[zone]
      tgt <- (wn$target_hr %||% 1) * cycle_hazard(s_c)[zone]
      above <- raw > tgt + 1e-12
      d$share_above_target <- mean(above)
      # Weighted by waning strength and by intervention patients still at risk:
      # a crossing in the far tail, where almost nobody is left, barely matters
      wt <- w[zone] * s_pre[zone]
      d$share_above_weighted <- if (sum(wt) > 0) sum(wt * above) / sum(wt) else 0
      d$first_above <- if (any(above)) tm[zone][which(above)[1]] else NA
      d$at_risk_first_above <- if (any(above)) s_pre[zone][which(above)[1]] else NA
    }
    out[[ep]] <- d
  }
  out
}

# Absolute OS gap (intervention - comparator) at landmark years
survival_gap <- function(curves, years, endpoint = "OS") {
  at <- \(s) stats::approx(curves$t, s, xout = years, rule = 2)$y
  data.frame(year = years, comp = at(curves$comp[[endpoint]]),
             int = at(curves$int[[endpoint]]),
             gap = at(curves$int[[endpoint]]) - at(curves$comp[[endpoint]]))
}

# Curves for every structure at the current timing (for overlay plots)
compare_structures <- function(pieces, spec, endpoint = "OS") {
  out <- lapply(waning_structures()$type, \(ty) {
    s2 <- spec
    s2$waning$type <- ty
    if (ty == "gradual" && !(is.finite(spec$waning$end) && spec$waning$end > spec$waning$start)) {
      s2$waning$end <- spec$waning$start + 2
    }
    cv <- assemble_curves(pieces, s2)
    data.frame(time = cv$t, surv = cv$int[[endpoint]], structure = waning_structures()$label[waning_structures()$type == ty])
  })
  cv0 <- assemble_curves(pieces, spec)
  rbind(do.call(rbind, out),
        data.frame(time = cv0$t, surv = cv0$comp[[endpoint]], structure = "Comparator"))
}

# Explorer: every structure x start time (TA798 / TA1099: test timings and structures)
waning_grid <- function(pieces, spec, starts = c(2, 3, 5, 7, 10), gradual_width = 2,
                        landmarks = c(5, 10)) {
  rows <- list()
  run <- function(ty, st) {
    s2 <- spec
    s2$waning$type <- ty
    s2$waning$start <- st
    s2$waning$end <- st + gradual_width
    cv <- assemble_curves(pieces, s2)
    r <- run_engine(cv, s2)
    g <- survival_gap(cv, landmarks)
    wt <- cv$waning
    data.frame(structure = waning_structures()$label[waning_structures()$type == ty],
               type = ty, start_input = if (ty == "none") NA else st,
               start_year = if (ty == "none") NA else wt$start,
               end_year = if (ty == "gradual") wt$end else if (ty == "none") NA else wt$start,
               inc_qaly = r$inc_qaly, inc_cost = r$inc_cost, icer = r$icer, nmb = r$nmb,
               t(setNames(g$int, paste0("int_os_yr", landmarks))),
               t(setNames(g$gap, paste0("gap_yr", landmarks))),
               check.names = FALSE)
  }
  rows[[1]] <- run("none", NA)
  for (ty in c("equal", "gradual", "stop_drop")) for (st in starts) rows[[length(rows) + 1]] <- run(ty, st)
  g <- do.call(rbind, rows)
  wtp <- spec$settings$wtp
  g$cost_effective <- (g$inc_qaly > 0 & is.finite(g$icer) & g$icer <= wtp) | (g$inc_qaly > 0 & g$inc_cost <= 0)
  g
}
