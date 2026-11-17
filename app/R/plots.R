# plots.R : ggplot2 helpers ------------------------------------------------

arm_colours <- function(arms) setNames(c("#3E5C8A", "#D9822B"), c(arms$comp, arms$int))

theme_psm <- function() {
  ggplot2::theme_minimal(base_size = 13) +
    ggplot2::theme(legend.position = "bottom", panel.grid.minor = ggplot2::element_blank())
}

pct_axis <- function() ggplot2::scale_y_continuous(labels = \(x) paste0(round(100 * x), "%"))

# KM + candidate fits for one arm, over the horizon
plot_km_fits <- function(fe, arm_key, models, horizon) {
  arm_label <- fe$arms[[arm_key]]
  km <- fe$km[[arm_key]]
  t <- seq(0, horizon, length.out = 400)
  fits <- fe$indep[[arm_key]][models]
  fits <- fits[vapply(fits, \(f) isTRUE(f$ok), logical(1))]
  curves <- do.call(rbind, lapply(fits, \(f) data.frame(time = t, surv = f$S(t, arm_label),
                                                        model = f$label)))
  ggplot2::ggplot() +
    ggplot2::geom_step(data = km_df(km, "KM"), ggplot2::aes(time, surv), linewidth = 0.9) +
    ggplot2::geom_line(data = curves, ggplot2::aes(time, surv, colour = model), linewidth = 0.7) +
    ggplot2::geom_vline(xintercept = max(km$time), linetype = 3) +
    pct_axis() + theme_psm() +
    ggplot2::labs(x = "Years", y = "Survival", colour = NULL,
                  title = paste(fe$endpoint, arm_label),
                  subtitle = "Black step: Kaplan-Meier. Dotted line: end of observed data.")
}

plot_fit_hazards <- function(fe, arm_key, models, horizon) {
  arm_label <- fe$arms[[arm_key]]
  t <- seq(0.05, horizon, length.out = 300)
  fits <- fe$indep[[arm_key]][models]
  fits <- fits[vapply(fits, \(f) isTRUE(f$ok), logical(1))]
  hz <- do.call(rbind, lapply(fits, \(f) data.frame(
    time = t, hazard = hazard_from_S(\(x) f$S(x, arm_label), t), model = f$label)))
  ggplot2::ggplot(hz, ggplot2::aes(time, hazard, colour = model)) +
    ggplot2::geom_line(linewidth = 0.7) +
    ggplot2::geom_vline(xintercept = max(fe$km[[arm_key]]$time), linetype = 3) +
    theme_psm() +
    ggplot2::labs(x = "Years", y = "Hazard (per year)", colour = NULL,
                  title = paste("Extrapolated hazards:", fe$endpoint, arm_label),
                  subtitle = "Choosing a distribution is choosing future mortality (Doc 7).")
}

plot_logcumhaz <- function(fe) {
  d <- do.call(rbind, lapply(names(fe$km), \(a) {
    km <- fe$km[[a]]
    keep <- km$surv > 0 & km$surv < 1 & km$time > 0
    data.frame(log_t = log(km$time[keep]), lch = log(-log(km$surv[keep])), arm = fe$arms[[a]])
  }))
  ggplot2::ggplot(d, ggplot2::aes(log_t, lch, colour = arm)) +
    ggplot2::geom_step(linewidth = 0.8) +
    ggplot2::scale_colour_manual(values = arm_colours(fe$arms)) + theme_psm() +
    ggplot2::labs(x = "log(time)", y = "log(-log S(t))", colour = NULL,
                  title = "Log-cumulative hazard",
                  subtitle = "Roughly parallel lines are compatible with proportional hazards.")
}

plot_schoenfeld <- function(zph) {
  d <- data.frame(time = zph$x, resid = zph$y[, 1])
  ggplot2::ggplot(d, ggplot2::aes(time, resid)) +
    ggplot2::geom_point(alpha = 0.35) +
    ggplot2::geom_smooth(method = "loess", formula = y ~ x, se = TRUE, colour = "#D9822B") +
    ggplot2::geom_hline(yintercept = 0, linetype = 2) + theme_psm() +
    ggplot2::labs(x = "Time (transformed)", y = "Scaled Schoenfeld residual, beta(t)",
                  title = "Scaled Schoenfeld residuals",
                  subtitle = "A trend over time suggests a time-varying treatment effect.")
}

# Final model curves (after effects, floor, constraints) with KM overlay
plot_model_curves <- function(curves, arms, fits = NULL, endpoints = c("OS", "PFS"),
                              show_gp = TRUE, xmax = NULL) {
  t <- curves$t
  d <- do.call(rbind, lapply(c("comp", "int"), \(a) do.call(rbind, lapply(endpoints, \(ep)
    data.frame(time = t, surv = curves[[a]][[ep]], arm = arms[[a]], endpoint = ep)))))
  p <- ggplot2::ggplot(d, ggplot2::aes(time, surv, colour = arm, linetype = endpoint)) +
    ggplot2::geom_line(linewidth = 0.8)
  if (!is.null(fits)) {
    kms <- do.call(rbind, lapply(endpoints, \(ep) {
      fe <- fits[[ep]]
      if (is.null(fe)) return(NULL)
      do.call(rbind, lapply(names(fe$km), \(a) {
        k <- km_df(fe$km[[a]], arms[[a]])
        k$endpoint <- ep
        k
      }))
    }))
    if (!is.null(kms)) {
      p <- p + ggplot2::geom_step(data = kms, ggplot2::aes(time, surv, colour = arm,
                                                          linetype = endpoint),
                                  alpha = 0.45, linewidth = 0.5)
    }
  }
  if (show_gp && !is.null(curves$gp)) {
    p <- p + ggplot2::geom_line(data = data.frame(time = t, surv = curves$gp),
                                ggplot2::aes(time, surv), inherit.aes = FALSE,
                                colour = "grey45", linetype = "dotdash")
  }
  p + ggplot2::scale_colour_manual(values = arm_colours(arms)) +
    ggplot2::coord_cartesian(xlim = c(0, xmax %||% max(t))) + pct_axis() + theme_psm() +
    ggplot2::labs(x = "Years", y = "Proportion of cohort", colour = NULL, linetype = NULL,
                  subtitle = "Thin steps: Kaplan-Meier. Grey dot-dash: general population.")
}

plot_model_hazards <- function(curves, arms, endpoint = "OS") {
  tm <- cycle_mid(curves$t)
  dt <- diff(curves$t)
  d <- do.call(rbind, lapply(c("comp", "int"), \(a)
    data.frame(time = tm, hazard = cycle_hazard(curves[[a]][[endpoint]]) / dt, arm = arms[[a]])))
  p <- ggplot2::ggplot(d, ggplot2::aes(time, hazard, colour = arm)) +
    ggplot2::geom_line(linewidth = 0.8)
  if (!is.null(curves$gp)) {
    p <- p + ggplot2::geom_line(data = data.frame(time = tm, hazard = cycle_hazard(curves$gp) / dt),
                                ggplot2::aes(time, hazard), inherit.aes = FALSE,
                                colour = "grey45", linetype = "dotdash")
  }
  p + ggplot2::scale_colour_manual(values = arm_colours(arms)) + theme_psm() +
    ggplot2::labs(x = "Years", y = "Hazard (per year)", colour = NULL,
                  title = paste("Modelled", endpoint, "hazard"),
                  subtitle = "Grey dot-dash: general-population hazard.")
}

plot_occupancy <- function(curves, arm_key, arm_label) {
  cv <- curves[[arm_key]]
  t <- curves$t
  d <- rbind(
    data.frame(time = t, prop = cv$PFS, state = "Progression-free"),
    data.frame(time = t, prop = cv$OS - cv$PFS, state = "Progressed"),
    data.frame(time = t, prop = 1 - cv$OS, state = "Dead")
  )
  d$state <- factor(d$state, levels = c("Dead", "Progressed", "Progression-free"))
  ggplot2::ggplot(d, ggplot2::aes(time, prop, fill = state)) +
    ggplot2::geom_area(alpha = 0.9) +
    ggplot2::geom_line(data = data.frame(time = t, prop = cv$ToT), ggplot2::aes(time, prop),
                       inherit.aes = FALSE, linetype = 2) +
    ggplot2::scale_fill_manual(values = c("Dead" = "#D3D1C7", "Progressed" = "#E8A87C",
                                          "Progression-free" = "#5DA9A4")) +
    pct_axis() + theme_psm() +
    ggplot2::labs(x = "Years", y = "Proportion of cohort", fill = NULL, title = arm_label,
                  subtitle = "Dashed line: proportion on treatment (ToT).")
}

plot_implied_hr <- function(curves, xmax = NULL) {
  d <- rbind(cbind(implied_hr(curves, "OS"), endpoint = "OS"),
             cbind(implied_hr(curves, "PFS"), endpoint = "PFS"))
  ggplot2::ggplot(d, ggplot2::aes(time, hr, colour = endpoint)) +
    ggplot2::geom_hline(yintercept = 1, linetype = 2) +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::coord_cartesian(xlim = c(0, xmax %||% max(d$time)),
                             ylim = c(0, max(1.5, stats::quantile(d$hr, 0.99, na.rm = TRUE)))) +
    theme_psm() +
    ggplot2::labs(x = "Years", y = "Hazard ratio, intervention vs comparator", colour = NULL,
                  title = "Implied hazard ratio over time",
                  subtitle = "Inspect the treatment effect the curves imply (Doc 3).")
}

plot_tornado <- function(dsa, outcome = c("icer", "nmb"), base_value) {
  outcome <- match.arg(outcome)
  lo <- dsa[[paste0(outcome, "_low")]]
  hi <- dsa[[paste0(outcome, "_high")]]
  d <- data.frame(parameter = dsa$parameter, lo = pmin(lo, hi), hi = pmax(lo, hi))
  d$range <- d$hi - d$lo
  d <- d[order(d$range), ]
  d$parameter <- factor(d$parameter, levels = d$parameter)
  ggplot2::ggplot(d) +
    ggplot2::geom_segment(ggplot2::aes(x = lo, xend = hi, y = parameter, yend = parameter),
                          linewidth = 6, colour = "#5DA9A4") +
    ggplot2::geom_vline(xintercept = base_value, linetype = 2) + theme_psm() +
    ggplot2::scale_x_continuous(labels = \(x) fmt_num(x, 0)) +
    ggplot2::labs(x = if (outcome == "icer") "ICER" else "Incremental net monetary benefit",
                  y = NULL, title = "One-way sensitivity analysis",
                  subtitle = "Dashed line: base case. Bars: result at low and high values.")
}

plot_ce_plane <- function(runs, wtp) {
  ok <- runs[runs$valid, ]
  ggplot2::ggplot(ok, ggplot2::aes(d_qaly, d_cost)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60") +
    ggplot2::geom_vline(xintercept = 0, colour = "grey60") +
    ggplot2::geom_point(alpha = 0.35, colour = "#3E5C8A") +
    ggplot2::geom_abline(slope = wtp, intercept = 0, linetype = 2) +
    ggplot2::geom_point(data = data.frame(d_qaly = mean(ok$d_qaly), d_cost = mean(ok$d_cost)),
                        colour = "#D9822B", size = 4) +
    ggplot2::scale_y_continuous(labels = \(x) fmt_num(x, 0)) + theme_psm() +
    ggplot2::labs(x = "Incremental QALYs", y = "Incremental cost",
                  title = "Cost-effectiveness plane",
                  subtitle = "Orange: mean of runs. Dashed: willingness-to-pay threshold.")
}

plot_ceac <- function(ce, wtp) {
  ggplot2::ggplot(ce, ggplot2::aes(wtp, p_ce)) +
    ggplot2::geom_line(linewidth = 0.9, colour = "#3E5C8A") +
    ggplot2::geom_vline(xintercept = wtp, linetype = 2) +
    ggplot2::scale_x_continuous(labels = \(x) fmt_num(x, 0)) +
    pct_axis() + theme_psm() +
    ggplot2::labs(x = "Willingness to pay per QALY", y = "Probability cost-effective",
                  title = "Cost-effectiveness acceptability curve")
}

# ---- V2: PSA survival uncertainty and convergence ---------------------------------

plot_psa_bands <- function(bands, curves, arms, xmax = 20) {
  det <- do.call(rbind, lapply(c("comp", "int"), \(a) do.call(rbind, lapply(c("OS", "PFS"), \(ep)
    data.frame(time = curves$t, surv = curves[[a]][[ep]], arm = arms[[a]], endpoint = ep)))))
  ggplot2::ggplot(bands, ggplot2::aes(time)) +
    ggplot2::geom_ribbon(ggplot2::aes(ymin = lower, ymax = upper, fill = arm), alpha = 0.25) +
    ggplot2::geom_line(ggplot2::aes(y = median, colour = arm), linewidth = 0.6, linetype = 2) +
    ggplot2::geom_line(data = det, ggplot2::aes(time, surv, colour = arm), linewidth = 0.9) +
    ggplot2::facet_wrap(~endpoint) +
    ggplot2::scale_colour_manual(values = arm_colours(arms)) +
    ggplot2::scale_fill_manual(values = arm_colours(arms)) +
    ggplot2::coord_cartesian(xlim = c(0, xmax)) + pct_axis() + theme_psm() +
    ggplot2::labs(x = "Years", y = "Proportion of cohort", colour = NULL, fill = NULL,
                  title = "Survival uncertainty across PSA runs",
                  subtitle = "Band: 95% of valid runs. Dashed: median run. Solid: deterministic curve.")
}

plot_convergence <- function(conv, outcome = c("icer", "nmb")) {
  outcome <- match.arg(outcome)
  d <- data.frame(runs = conv$runs, value = conv[[outcome]])
  final <- utils::tail(d$value, 1)
  ggplot2::ggplot(d, ggplot2::aes(runs, value)) +
    ggplot2::geom_hline(yintercept = final, linetype = 2, colour = "grey50") +
    ggplot2::geom_line(linewidth = 0.8, colour = "#3E5C8A") +
    ggplot2::scale_y_continuous(labels = \(x) fmt_num(x, 0)) + theme_psm() +
    ggplot2::labs(x = "Number of valid runs", y = if (outcome == "icer") "Running ICER (ratio of means)" else "Running net monetary benefit",
                  title = "Convergence",
                  subtitle = "A flat line by the end suggests enough runs. Dashed: final value.")
}

# ---- V3: waning ---------------------------------------------------------------

structure_colours <- function() {
  c("Comparator" = "#2C2C2A", "Sustained effect" = "#3E5C8A", "Equal hazards from start" = "#1D9E75",
    "Gradual equalisation" = "#D9822B", "Stop and drop" = "#A32D2D")
}

plot_waning_timeline <- function(cap, fu, start, end, type, xmax) {
  rows <- list()
  add <- \(lane, from, to, what) if (is.finite(from) && is.finite(to) && to > from) {
    rows[[length(rows) + 1]] <<- data.frame(lane = lane, from = from, to = min(to, xmax), what = what)
  }
  add("Treatment", 0, cap, "On treatment")
  add("Evidence", 0, fu, "Observed follow-up")
  add("Evidence", fu, xmax, "Extrapolated")
  if (type == "none") {
    add("Effect", 0, xmax, "Full effect")
  } else {
    add("Effect", 0, start, "Full effect")
    if (type == "gradual") add("Effect", start, end, "Waning")
    add("Effect", if (type == "gradual") end else start, xmax,
        if (type == "stop_drop") "Gain removed" else "Effect waned (gain kept)")
  }
  d <- do.call(rbind, rows)
  d$lane <- factor(d$lane, levels = c("Effect", "Evidence", "Treatment"))
  fill <- c("On treatment" = "#5DA9A4", "Observed follow-up" = "#3E5C8A", "Extrapolated" = "#D3D1C7",
            "Full effect" = "#1D9E75", "Waning" = "#EF9F27", "Effect waned (gain kept)" = "#F5C4B3",
            "Gain removed" = "#A32D2D")
  ggplot2::ggplot(d) +
    ggplot2::geom_rect(ggplot2::aes(xmin = from, xmax = to, ymin = as.numeric(lane) - 0.35,
                                    ymax = as.numeric(lane) + 0.35, fill = what)) +
    ggplot2::scale_y_continuous(breaks = 1:3, labels = levels(d$lane)) +
    ggplot2::scale_fill_manual(values = fill) +
    ggplot2::coord_cartesian(xlim = c(0, xmax)) + theme_psm() +
    ggplot2::theme(panel.grid.major.y = ggplot2::element_blank()) +
    ggplot2::labs(x = "Years since treatment start", y = NULL, fill = NULL)
}

plot_structures_compare <- function(df, current, xmax = 15, endpoint = "OS") {
  df$current <- df$structure == current | df$structure == "Comparator"
  ggplot2::ggplot(df, ggplot2::aes(time, surv, colour = structure)) +
    ggplot2::geom_line(ggplot2::aes(linewidth = current)) +
    ggplot2::scale_linewidth_manual(values = c(`TRUE` = 1.3, `FALSE` = 0.6), guide = "none") +
    ggplot2::scale_colour_manual(values = structure_colours()) +
    ggplot2::coord_cartesian(xlim = c(0, xmax)) + pct_axis() + theme_psm() +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2)) +
    ggplot2::labs(x = "Years since treatment start", y = paste("Intervention", endpoint), colour = NULL,
                  subtitle = "Thick line: current choice. Same timing for every structure.")
}

plot_gap_compare <- function(df, xmax = 15, landmarks = c(3, 5, 10)) {
  comp <- df[df$structure == "Comparator", ]
  g <- df[df$structure != "Comparator", ]
  g$gap <- g$surv - comp$surv[match(g$time, comp$time)]
  pts <- do.call(rbind, lapply(split(g, g$structure), \(x) {
    data.frame(time = landmarks, gap = stats::approx(x$time, x$gap, xout = landmarks)$y,
               structure = x$structure[1])
  }))
  ggplot2::ggplot(g, ggplot2::aes(time, gap, colour = structure)) +
    ggplot2::geom_hline(yintercept = 0, colour = "grey60") +
    ggplot2::geom_line(linewidth = 0.8) +
    ggplot2::geom_point(data = pts, size = 2) +
    ggplot2::geom_text(data = pts, ggplot2::aes(label = round(100 * gap, 1)), vjust = -0.8, size = 3.2,
                       show.legend = FALSE, check_overlap = TRUE) +
    ggplot2::scale_colour_manual(values = structure_colours()) +
    ggplot2::scale_y_continuous(labels = \(x) paste(round(100 * x, 1), "pts")) +
    ggplot2::coord_cartesian(xlim = c(0, xmax)) + theme_psm() +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2)) +
    ggplot2::labs(x = "Years since treatment start", y = "Absolute survival advantage", colour = NULL,
                  subtitle = "Equal hazards are not equal survival: the gain narrows.")
}

plot_waning_explorer <- function(g, wtp, outcome = c("icer", "inc_qaly", "nmb"), currency = "") {
  outcome <- match.arg(outcome)
  lab <- c(icer = paste("ICER", currency), inc_qaly = "Incremental QALYs", nmb = paste("Net monetary benefit", currency))
  base <- g[g$type == "none", ]
  d <- g[g$type != "none", ]
  p <- ggplot2::ggplot(d, ggplot2::aes(start_year, .data[[outcome]], colour = structure)) +
    ggplot2::geom_hline(yintercept = base[[outcome]], colour = structure_colours()[["Sustained effect"]],
                        linetype = 2) +
    ggplot2::geom_line(linewidth = 0.8) + ggplot2::geom_point(size = 2)
  if (outcome == "icer") p <- p + ggplot2::geom_hline(yintercept = wtp, colour = "grey40", linetype = 3)
  if (outcome == "nmb") p <- p + ggplot2::geom_hline(yintercept = 0, colour = "grey40", linetype = 3)
  p + ggplot2::scale_colour_manual(values = structure_colours()) +
    ggplot2::scale_y_continuous(labels = \(x) fmt_num(x, if (outcome == "inc_qaly") 2 else 0)) + theme_psm() +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 1)) +
    ggplot2::labs(x = "Year waning starts (since treatment start)", y = lab[[outcome]], colour = NULL,
                  subtitle = paste("Dashed: sustained effect.", if (outcome == "icer") "Dotted: threshold." else if (outcome == "nmb") "Dotted: break-even." else ""))
}
