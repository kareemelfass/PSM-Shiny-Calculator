# fct_survival.R : fitting, prediction, KM, hybrid curves ----------------------
# One adapter (fit_one) hides the fitting engine. Everything downstream only
# uses fit$S(t, arm), so the rest of the app is engine-agnostic.

model_catalogue <- function() {
  data.frame(
    id = c("exp", "weibull", "gompertz", "lnorm", "llogis", "gamma", "gengamma",
           "spline1", "spline2", "spline3"),
    label = c("Exponential", "Weibull", "Gompertz", "Log-normal", "Log-logistic",
              "Gamma", "Generalised gamma", "RP spline (1 knot)",
              "RP spline (2 knots)", "RP spline (3 knots)"),
    survreg = c("exponential", "weibull", NA, "lognormal", "loglogistic",
                NA, NA, NA, NA, NA)
  )
}

fit_engine <- function() {
  if (requireNamespace("flexsurv", quietly = TRUE)) "flexsurv" else "survreg"
}

available_models <- function(engine = fit_engine()) {
  cat <- model_catalogue()
  if (engine == "flexsurv") cat else cat[!is.na(cat$survreg), ]
}

model_label <- function(id) {
  cat <- model_catalogue()
  cat$label[match(id, cat$id)]
}

# ---- engine adapters --------------------------------------------------------

fit_flexsurv <- function(d, model, covariate) {
  f <- if (covariate) survival::Surv(time, event) ~ arm else survival::Surv(time, event) ~ 1
  fit <- if (startsWith(model, "spline")) {
    flexsurv::flexsurvspline(f, data = d, k = as.integer(sub("spline", "", model)),
                             scale = "hazard")
  } else {
    flexsurv::flexsurvreg(f, data = d, dist = model)
  }
  lev <- levels(d$arm)
  S <- function(t, arm = NULL) {
    out <- rep(1, length(t))
    pos <- t > 0
    if (any(pos)) {
      nd <- if (covariate) data.frame(arm = factor(arm, levels = lev)) else NULL
      s <- summary(fit, newdata = nd, t = t[pos], type = "survival", ci = FALSE, tidy = TRUE)
      est <- if ("est" %in% names(s)) s$est else s[[2]]
      out[pos] <- est[match(t[pos], s$time)]
    }
    out
  }

  # Survival from one set of natural-scale parameters (covariate effects already applied)
  S_nat <- function(par, t) {
    out <- rep(1, length(t))
    pos <- t > 0
    if (any(pos)) {
      args <- c(list(q = t[pos], lower.tail = FALSE), as.list(par), fit$aux)
      out[pos] <- do.call(fit$dfns$p, args)
    }
    out
  }

  # PSA: n draws from the asymptotic normal distribution of the parameters.
  # One draw serves both arms of a joint model, so their correlation is kept.
  draw_mvn <- function(n) {
    nd <- if (covariate) data.frame(arm = factor(lev, levels = lev)) else NULL
    sims <- flexsurv::normboot.flexsurvreg(fit, B = n, newdata = nd)
    if (!is.list(sims)) sims <- list(sims)
    names(sims) <- if (covariate) lev else "all"

    # Self-check: the parameter route must reproduce the fitted curve
    ref_par <- fit$res[fit$dlist$pars, "est"]
    tt <- c(0.5, 1, 2, 5)
    gap <- max(abs(S_nat(ref_par, tt) - S(tt, lev[1])))
    if (!is.finite(gap) || gap > 1e-6) {
      stop("flexsurv parameter sampling does not reproduce the fitted ", model,
           " curve (difference ", signif(gap, 3), "). Use the bootstrap method instead.")
    }
    function(i, t, arm = NULL) {
      m <- if (covariate) sims[[as.character(arm)]] else sims[[1]]
      S_nat(m[i, ], t)
    }
  }
  list(npar = fit$npars, loglik = fit$loglik, S = S, draw_mvn = draw_mvn, object = fit)
}

# n draws from a multivariate normal (base R; eigen fallback if not positive definite)
mvn_draw <- function(n, mu, V) {
  k <- length(mu)
  R <- tryCatch(chol(V), error = \(e) {
    ev <- eigen(V, symmetric = TRUE)
    t(ev$vectors %*% diag(sqrt(pmax(ev$values, 0)), k))
  })
  z <- matrix(stats::rnorm(n * k), n, k)
  sweep(z %*% R, 2, mu, "+")
}

fit_survreg <- function(d, model, covariate) {
  dist <- model_catalogue()$survreg[model_catalogue()$id == model]
  if (is.na(dist)) stop(model, " needs the flexsurv package")
  f <- if (covariate) survival::Surv(time, event) ~ arm else survival::Surv(time, event) ~ 1
  fit <- survival::survreg(f, data = d, dist = dist)
  b <- stats::coef(fit)
  nb <- length(b)
  ref <- levels(d$arm)[1]

  # theta = regression coefficients, then log(scale) unless exponential.
  # This matches the order and scale of vcov(fit).
  theta <- c(b, if (dist != "exponential") log(fit$scale))
  V <- stats::vcov(fit)
  if (nrow(V) != length(theta)) stop("Unexpected survreg covariance dimensions")

  S_par <- function(th, t, arm = NULL) {
    sc <- if (dist == "exponential") 1 else exp(th[[nb + 1]])
    lp <- th[[1]] + if (covariate) th[[2]] * as.numeric(arm != ref) else 0
    z <- (log(pmax(t, 1e-12)) - lp) / sc
    s <- switch(dist,
      exponential = , weibull = exp(-exp(z)),
      lognormal   = stats::pnorm(z, lower.tail = FALSE),
      loglogistic = stats::plogis(z, lower.tail = FALSE)
    )
    s[t <= 0] <- 1
    s
  }
  S <- function(t, arm = NULL) S_par(theta, t, arm)

  # PSA: n draws of theta from N(theta_hat, vcov). One draw serves both arms.
  draw_mvn <- function(n) {
    th <- mvn_draw(n, theta, V)
    function(i, t, arm = NULL) S_par(th[i, ], t, arm)
  }
  list(npar = length(theta), loglik = fit$loglik[2], S = S, draw_mvn = draw_mvn,
       theta = theta, vcov = V, object = fit)
}

# Fit one model. d: data.frame(time, event, arm factor with comparator first)
fit_one <- function(d, model, covariate = FALSE, engine = fit_engine()) {
  res <- tryCatch(
    if (engine == "flexsurv") fit_flexsurv(d, model, covariate) else fit_survreg(d, model, covariate),
    error = \(e) list(error = conditionMessage(e))
  )
  if (!is.null(res$error) || !is.finite(res$loglik %||% NA)) {
    return(list(ok = FALSE, model = model, error = res$error %||% "did not converge"))
  }
  n <- nrow(d)
  c(res, list(
    ok = TRUE, model = model, label = model_label(model), engine = engine,
    n = n, events = sum(d$event),
    aic = -2 * res$loglik + 2 * res$npar,
    bic = -2 * res$loglik + res$npar * log(n)   # BIC with n = patients
  ))
}

# Fit all models for one endpoint: independent per arm + joint (arm covariate)
fit_endpoint <- function(ipd, endpoint, arms, engine = fit_engine(), progress = NULL) {
  d_all <- ipd[ipd$endpoint == endpoint, c("arm", "time", "event")]
  d_all <- d_all[d_all$arm %in% unname(unlist(arms)), ]
  d_all$arm <- factor(d_all$arm, levels = c(arms$comp, arms$int))
  models <- available_models(engine)$id
  have <- c(comp = any(d_all$arm == arms$comp), int = any(d_all$arm == arms$int))

  out <- list(endpoint = endpoint, arms = arms, engine = engine, data = d_all,
              have = have, indep = list(), joint = list(), km = list())
  for (a in names(have)[have]) {
    d <- d_all[d_all$arm == arms[[a]], ]
    d$arm <- droplevels(d$arm)
    out$km[[a]] <- survival::survfit(survival::Surv(time, event) ~ 1, data = d)
    for (m in models) {
      if (!is.null(progress)) progress(paste(endpoint, arms[[a]], model_label(m)))
      out$indep[[a]][[m]] <- fit_one(d, m, FALSE, engine)
    }
  }
  if (all(have)) {
    for (m in models) {
      if (!is.null(progress)) progress(paste(endpoint, "joint", model_label(m)))
      out$joint[[m]] <- fit_one(d_all, m, TRUE, engine)
    }
  }
  out$gof <- gof_table(out)
  out
}

gof_table <- function(fe) {
  rows <- list()
  add <- function(approach, arm, f) {
    rows[[length(rows) + 1]] <<- data.frame(
      approach = approach, arm = arm, model = f$model, label = model_label(f$model),
      npar = if (f$ok) f$npar else NA, loglik = if (f$ok) f$loglik else NA,
      aic = if (f$ok) f$aic else NA, bic = if (f$ok) f$bic else NA,
      note = if (f$ok) "" else f$error
    )
  }
  for (a in names(fe$indep)) for (f in fe$indep[[a]]) add("Independent", fe$arms[[a]], f)
  for (f in fe$joint) add("Joint (arm covariate)", "Both arms", f)
  if (!length(rows)) return(data.frame())
  g <- do.call(rbind, rows)
  grp <- paste(g$approach, g$arm)
  g$d_aic <- g$aic - stats::ave(g$aic, grp, FUN = \(x) min(x, na.rm = TRUE))
  g$d_bic <- g$bic - stats::ave(g$bic, grp, FUN = \(x) min(x, na.rm = TRUE))
  g
}

# Kaplan-Meier survival at times t (step function, carried forward)
km_surv <- function(km, t) {
  s <- summary(km, times = t, extend = TRUE)$surv
  s[t <= 0] <- 1
  s
}

km_at_risk <- function(km, t) summary(km, times = t, extend = TRUE)$n.risk

km_df <- function(km, arm_label) {
  data.frame(time = c(0, km$time), surv = c(1, km$surv), arm = arm_label)
}

# Survival from a fit, optionally hybrid: KM up to cut, then conditional fit
fit_curve <- function(fit, arm_label, t, km = NULL, hybrid = FALSE, cut = NA) {
  hybrid_curve(\(x) fit$S(x, arm_label), t, km, hybrid, cut)
}

# Same rule for any survival function (used by the PSA with sampled curves)
hybrid_curve <- function(S_fun, t, km = NULL, hybrid = FALSE, cut = NA) {
  s <- S_fun(t)
  if (isTRUE(hybrid) && !is.null(km) && is.finite(cut)) {
    s_km <- km_surv(km, t)
    s_cut <- km_surv(km, cut)
    s <- ifelse(t <= cut, s_km, s_cut * s / S_fun(cut))
  }
  s
}

# Numeric hazard (per year) from a survival function, for plots
hazard_from_S <- function(S_fun, t, eps = 1e-3) {
  t <- pmax(t, eps)
  (log(pmax(S_fun(t), 1e-300)) - log(pmax(S_fun(t + eps), 1e-300))) / eps
}

# Cox HR (intervention vs comparator) and PH test for one endpoint
cox_summary <- function(fe) {
  if (!all(fe$have)) return(NULL)
  d <- fe$data
  fit <- survival::coxph(survival::Surv(time, event) ~ arm, data = d)
  zph <- tryCatch(survival::cox.zph(fit), error = \(e) NULL)
  list(
    hr = exp(stats::coef(fit))[[1]],
    se_log = sqrt(diag(stats::vcov(fit)))[[1]],
    p_ph = if (is.null(zph)) NA else zph$table[1, "p"],
    zph = zph
  )
}
