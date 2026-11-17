# fct_effects.R : from fits to final model curves -------------------------------
# Stage 1 (predict_pieces): evaluate fitted curves on the grid. Slow, run once.
# Stage 2 (assemble_curves): HR, waning, gen-pop floor, constraints. Fast, so
# DSA and PSA can rerun it many times with different HRs.

# ---- hazard ratios ----------------------------------------------------------

# Convert a reported HR to "intervention vs comparator" on the log scale.
# direction = "int_vs_comp" (HR < 1 favours intervention) or "comp_vs_int".
hr_to_log <- function(est, lower, upper, direction = "int_vs_comp") {
  if (direction == "comp_vs_int") {
    est <- 1 / est
    tmp <- 1 / upper
    upper <- 1 / lower
    lower <- tmp
  }
  se <- if (is.finite(lower) && is.finite(upper) && lower > 0 && upper > lower) {
    (log(upper) - log(lower)) / (2 * stats::qnorm(0.975))
  } else NA_real_
  list(log_hr = log(est), se = se, hr = est, lower = lower, upper = upper)
}

# ---- waning (V3) ------------------------------------------------------------
# Structures (Doc 12):
#   none      sustained effect: the fitted/HR-based intervention curve is kept
#   equal     hazards converge to target_hr x comparator hazard at `start`
#             (immediate hazard equalisation; survival already gained is kept)
#   gradual   the relative effect moves linearly to target_hr between `start` and `end`
#   stop_drop intervention survival set equal to the comparator's from `start`
#             (an instantaneous loss of the gain: a stress test, not a finite hazard)
# target_hr  hazard ratio reached after full waning (1 = no relative effect left)
# frailty    gamma frailty variance. 0 = wane the population (marginal) hazards, as in V2.
#            > 0 = wane the individual (conditional) hazards: the arms' surviving
#            populations differ in frailty, so after waning the marginal HR goes above
#            the target (Jennings et al. 2024).

# Weight w(t): 0 = full intervention effect, 1 = effect fully waned
waning_weight <- function(tm, start, end) {
  if (!is.finite(end) || end <= start) return(as.numeric(tm >= start))
  pmin(pmax((tm - start) / (end - start), 0), 1)
}

# Gamma frailty: marginal <-> conditional cumulative hazard (identity when theta = 0)
to_conditional <- function(H, theta) if (theta > 0) expm1(pmin(theta * H, 700)) / theta else H
to_marginal <- function(H, theta) if (theta > 0) log1p(theta * H) / theta else H

apply_waning <- function(s_int, s_comp, t, type = "none", start = Inf, end = Inf,
                         target_hr = 1, frailty = 0) {
  if (type == "none" || !is.finite(start)) return(s_int)
  if (type == "stop_drop") {
    out <- s_int
    out[t >= start] <- s_comp[t >= start]
    return(out)
  }
  if (type == "equal") end <- start
  theta <- max(frailty %||% 0, 0)
  cum_h <- \(s) -log(pmax(s, 1e-300))
  dh_i <- diff(to_conditional(cum_h(s_int), theta))   # per-cycle hazards on the chosen scale
  dh_c <- diff(to_conditional(cum_h(s_comp), theta))
  w <- waning_weight(cycle_mid(t), start, end)
  H <- c(0, cumsum((1 - w) * dh_i + w * target_hr * dh_c))
  exp(-to_marginal(H, theta))
}

# Absolute waning times. anchor "start": years since treatment start (model start).
# anchor "tx_stop": years after the intervention's maximum treatment duration.
waning_times <- function(spec) {
  wn <- spec$waning
  cap <- spec$tot$int$cap %||% NA
  use_cap <- identical(wn$anchor, "tx_stop") && is.finite(cap) && cap > 0
  shift <- if (use_cap) cap else 0
  list(start = wn$start + shift,
       end = if (wn$type == "gradual") wn$end + shift else wn$start + shift,
       anchored_to_cap = use_cap,
       anchor_missing_cap = identical(wn$anchor, "tx_stop") && !use_cap)
}

# ---- general-population mortality --------------------------------------------

# Cumulative hazard from the life table's first age to age a (piecewise constant)
lt_cum_hazard <- function(lt_age, qx, a) {
  mu <- -log(1 - pmin(qx, 0.999999))
  width <- c(diff(lt_age), 1)
  h_start <- c(0, cumsum(mu * width))[seq_along(lt_age)]
  i <- pmax(findInterval(a, lt_age), 1)
  h_start[i] + mu[i] * (a - lt_age[i])
}

# General-population survival on the grid, mixed by sex
gen_pop_surv <- function(t, lifetable, age0, male_share) {
  lt <- lifetable[order(lifetable$age), ]
  s_sex <- function(q) {
    exp(-(lt_cum_hazard(lt$age, q, age0 + t) - lt_cum_hazard(lt$age, q, age0)))
  }
  male_share * s_sex(lt$qx_male) + (1 - male_share) * s_sex(lt$qx_female)
}

# ---- stage 1: predict pieces ------------------------------------------------

# spec$endpoints$OS / $PFS: approach ("indep"|"joint"|"hr"), comp_model,
#   int_model, hybrid, cut
# spec$tot$comp / $int: method ("fit"|"pfs"), model
predict_pieces <- function(spec, fits, t) {
  arms <- spec$arms
  one_endpoint <- function(ep) {
    es <- spec$endpoints[[ep]]
    fe <- fits[[ep]]
    get_fit <- function(a, model) {
      f <- if (es$approach == "joint") fe$joint[[model]] else fe$indep[[a]][[model]]
      if (is.null(f) || !isTRUE(f$ok)) stop(ep, ": model '", model_label(model),
                                            "' is not available for ", arms[[a]])
      f
    }
    comp <- fit_curve(get_fit("comp", es$comp_model), arms$comp, t,
                      fe$km$comp, es$hybrid, es$cut)
    int <- switch(es$approach,
      hr = NULL,
      joint = fit_curve(get_fit("int", es$comp_model), arms$int, t,
                        fe$km$int, es$hybrid, es$cut),
      indep = fit_curve(get_fit("int", es$int_model), arms$int, t,
                        fe$km$int, es$hybrid, es$cut)
    )
    list(comp = comp, int = int)
  }
  tot_piece <- function(a) {
    ts <- spec$tot[[a]]
    if (ts$method != "fit") return(NULL)
    f <- fits$ToT$indep[[a]][[ts$model]]
    if (is.null(f) || !isTRUE(f$ok)) stop("ToT model not available for ", arms[[a]])
    f$S(t, arms[[a]])
  }
  st <- spec$settings
  list(
    t = t,
    OS = one_endpoint("OS"),
    PFS = one_endpoint("PFS"),
    ToT = list(comp = tot_piece("comp"), int = tot_piece("int")),
    gp = if (isTRUE(st$floor_on) && !is.null(st$lifetable)) {
      gen_pop_surv(t, st$lifetable, st$age, st$male)
    } else NULL
  )
}

# ---- stage 2: assemble final curves -------------------------------------------

# hr: named list(OS = , PFS = ) of HRs (int vs comp), used where approach == "hr"
assemble_curves <- function(pieces, spec, hr = NULL) {
  t <- pieces$t
  st <- spec$settings
  wn <- spec$waning
  wt <- waning_times(spec)
  hr <- hr %||% list(OS = spec$endpoints$OS$hr$hr, PFS = spec$endpoints$PFS$hr$hr)
  flags <- list()

  curves <- list(comp = list(), int = list())
  pre_waning <- list()
  for (ep in c("OS", "PFS")) {
    p <- pieces[[ep]]
    s_comp <- p$comp
    s_int <- if (spec$endpoints[[ep]]$approach == "hr") s_comp^hr[[ep]] else p$int
    pre_waning[[ep]] <- s_int
    apply_here <- isTRUE(wn[[paste0("apply_", tolower(ep))]])
    if (apply_here) {
      s_int <- apply_waning(s_int, s_comp, t, wn$type, wt$start, wt$end,
                            wn$target_hr %||% 1, wn$frailty %||% 0)
    }
    curves$comp[[ep]] <- s_comp
    curves$int[[ep]] <- s_int
  }

  # General-population mortality floor on OS and PFS hazards
  gp_dh <- if (!is.null(pieces$gp)) cycle_hazard(pieces$gp) else NULL
  for (a in c("comp", "int")) {
    flags[[a]] <- list()
    if (!is.null(gp_dh)) {
      dh_os <- cycle_hazard(curves[[a]]$OS)
      binds <- which(gp_dh > dh_os)
      flags[[a]]$floor_from <- if (length(binds)) t[min(binds)] else NA
      flags[[a]]$floor_cycles <- length(binds)
      curves[[a]]$OS <- from_hazard(pmax(dh_os, gp_dh))
      curves[[a]]$PFS <- from_hazard(pmax(cycle_hazard(curves[[a]]$PFS), gp_dh))
    }
  }

  # Explicit, reported constraints (never silent)
  tol <- 1e-10
  for (a in c("comp", "int")) {
    cv <- curves[[a]]
    cross <- which(cv$PFS > cv$OS + tol)
    flags[[a]]$pfs_gt_os_cycles <- length(cross)
    flags[[a]]$pfs_gt_os_from <- if (length(cross)) t[min(cross)] else NA
    if (isTRUE(st$pfs_le_os)) cv$PFS <- pmin(cv$PFS, cv$OS)

    ts <- spec$tot[[a]]
    tot <- if (ts$method == "fit") pieces$ToT[[a]] else cv$PFS
    if (is.finite(ts$cap %||% NA) && ts$cap > 0) tot[t > ts$cap] <- 0
    if (!isTRUE(ts$beyond_prog)) {
      over_pfs <- which(tot > cv$PFS + tol)
      flags[[a]]$tot_gt_pfs_cycles <- length(over_pfs)
      tot <- pmin(tot, cv$PFS)
    }
    over_os <- which(tot > cv$OS + tol)
    flags[[a]]$tot_gt_os_cycles <- length(over_os)
    flags[[a]]$tot_gt_os_from <- if (length(over_os)) t[min(over_os)] else NA
    cv$ToT <- pmin(tot, cv$OS)
    curves[[a]] <- cv
  }

  list(t = t, comp = curves$comp, int = curves$int, gp = pieces$gp,
       hr_used = hr, flags = flags, pre_waning = pre_waning,
       waning = utils::modifyList(wn, wt))   # absolute start/end replace the user offsets
}

# Implied hazard ratio per cycle (intervention / comparator)
implied_hr <- function(curves, endpoint) {
  hi <- cycle_hazard(curves$int[[endpoint]])
  hc <- cycle_hazard(curves$comp[[endpoint]])
  data.frame(time = cycle_mid(curves$t), hr = ifelse(hc > 1e-12, hi / hc, NA))
}
