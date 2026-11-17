# fct_engine.R : state occupancy -> life-years, QALYs, costs --------------------
# Every formula is here, in plain R, so each number in the app can be traced.

# One arm. cv = list(OS, PFS, ToT) on grid t; e = arm economics; s = settings
arm_outcomes <- function(t, cv, e, u, s) {
  e <- lapply(e, \(x) unname(as.numeric(x)))
  u <- lapply(u, \(x) unname(as.numeric(x)))
  dt <- diff(t)
  tm <- cycle_mid(t)
  pf <- cv$PFS
  pd <- cv$OS - cv$PFS
  dead <- 1 - cv$OS

  # Time in a state during each cycle: trapezoid (half-cycle correction) or start
  in_cycle <- function(x) {
    if (isTRUE(s$hcc)) dt * (x[-length(x)] + x[-1]) / 2 else dt * x[-length(x)]
  }
  df_out <- 1 / (1 + s$dr_out)^tm
  df_cost <- 1 / (1 + s$dr_cost)^tm

  ly_pf <- in_cycle(pf)
  ly_pd <- in_cycle(pd)
  on_tx <- in_cycle(cv$ToT)                     # expected years on treatment
  new_prog <- -diff(cv$PFS) * e$prog_share      # entries to progressed disease
  new_death <- -diff(cv$OS)

  per_year_drug <- e$drug_cost * (365.25 / e$interval_days) * e$rdi
  per_year_admin <- e$admin_cost * (365.25 / e$interval_days)

  cost_cycle <- data.frame(
    drug = on_tx * per_year_drug,
    admin = on_tx * per_year_admin,
    pf_state = ly_pf * e$c_pf,
    pd_state = ly_pd * e$c_pd,
    subsequent = new_prog * e$c_subseq,
    terminal = new_death * e$c_terminal
  )
  qaly_cycle <- u$pf * ly_pf + u$pd * ly_pd

  trace <- data.frame(
    cycle = seq_along(tm), t_start = t[-length(t)], t_end = t[-1],
    pf = pf[-1], pd = pd[-1], dead = dead[-1], on_tx_prop = cv$ToT[-1],
    ly = ly_pf + ly_pd, qaly = qaly_cycle,
    cost_cycle, disc_out = df_out, disc_cost = df_cost
  )

  disc_costs <- colSums(cost_cycle * df_cost)
  one_off <- c(upfront = e$upfront, ae = e$ae_cost)
  costs <- c(one_off["upfront"], disc_costs, one_off["ae"])

  list(
    trace = trace,
    summary = c(
      ly_undisc = sum(ly_pf + ly_pd),
      ly_disc = sum((ly_pf + ly_pd) * df_out),
      ly_pf = sum(ly_pf * df_out),
      ly_pd = sum(ly_pd * df_out),
      years_on_tx = sum(on_tx),
      qaly_undisc = sum(qaly_cycle) - e$ae_qaly,
      qaly = sum(qaly_cycle * df_out) - e$ae_qaly,
      setNames(costs, paste0("cost_", names(costs))),
      cost_total = sum(costs)
    )
  )
}

icer_label <- function(d_cost, d_qaly) {
  if (!is.finite(d_qaly) || abs(d_qaly) < 1e-12) return("Undefined (no QALY difference)")
  if (d_qaly > 0 && d_cost <= 0) return("Dominant (more effective, not more costly)")
  if (d_qaly < 0 && d_cost >= 0) return("Dominated (less effective, not less costly)")
  if (d_qaly < 0) return("South-west quadrant: cost saved per QALY lost")
  "Cost per QALY gained"
}

# Full model run
run_engine <- function(curves, spec) {
  s <- spec$settings
  u <- spec$econ$utility
  res <- lapply(c(comp = "comp", int = "int"),
                \(a) arm_outcomes(curves$t, curves[[a]], spec$econ[[a]], u, s))
  comp <- res$comp$summary
  int <- res$int$summary
  inc <- int - comp

  # Share of incremental (undiscounted) life-years accrued after trial follow-up
  d_ly <- res$int$trace$ly - res$comp$trace$ly
  fu <- s$fu_years %||% NA
  share_after <- if (is.finite(fu) && sum(d_ly) > 0) {
    sum(d_ly[res$int$trace$t_start >= fu]) / sum(d_ly)
  } else NA

  icer <- if (abs(inc[["qaly"]]) > 1e-12) inc[["cost_total"]] / inc[["qaly"]] else NA
  list(
    comp = res$comp, int = res$int,
    table = data.frame(measure = names(comp), comparator = comp,
                       intervention = int, incremental = inc, row.names = NULL),
    inc_cost = inc[["cost_total"]], inc_qaly = inc[["qaly"]], inc_ly = inc[["ly_disc"]],
    icer = icer, icer_label = icer_label(inc[["cost_total"]], inc[["qaly"]]),
    nmb = s$wtp * inc[["qaly"]] - inc[["cost_total"]],
    share_ly_after_fu = share_after
  )
}

# Pretty results table for display
results_display <- function(r, arms, currency) {
  lab <- c(
    ly_undisc = "Life-years (undiscounted)", ly_disc = "Life-years (discounted)",
    ly_pf = "\u00a0\u00a0\u00a0of which progression-free", ly_pd = "\u00a0\u00a0\u00a0of which progressed",
    years_on_tx = "Years on treatment (undiscounted)",
    qaly_undisc = "QALYs (undiscounted)", qaly = "QALYs (discounted)",
    cost_upfront = "Upfront cost", cost_drug = "Drug acquisition",
    cost_admin = "Administration", cost_pf_state = "Progression-free care",
    cost_pd_state = "Progressed-disease care", cost_subsequent = "Subsequent therapy",
    cost_terminal = "Terminal care", cost_ae = "Adverse events",
    cost_total = "Total cost (discounted)"
  )
  tb <- r$table
  is_cost <- startsWith(tb$measure, "cost_")
  fmt <- \(x) ifelse(is_cost, fmt_num(x, 0), fmt_num(x, 3))
  out <- data.frame(
    Measure = lab[tb$measure],
    a = fmt(tb$comparator), b = fmt(tb$intervention), c = fmt(tb$incremental)
  )
  names(out) <- c("Measure", arms$comp, arms$int, "Incremental")
  out$Measure[is_cost] <- paste0(out$Measure[is_cost], " (", currency, ")")
  out
}
