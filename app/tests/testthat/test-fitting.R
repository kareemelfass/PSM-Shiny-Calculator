# Fitting, full pipeline on the example data, and the flexsurv cross-check

test_that("models fit on the example data and the GOF table is complete", {
  fe <- fit_endpoint(example_ipd, "OS", example_arms)
  n_models <- nrow(available_models())
  expect_equal(nrow(fe$gof), 3 * n_models)            # 2 arms + joint
  expect_true(all(fe$gof$d_aic >= 0, na.rm = TRUE))
  w <- fe$indep$comp$weibull
  expect_true(w$ok)
  expect_equal(w$S(0), 1)
  expect_true(all(diff(w$S(seq(0, 10, 0.5))) <= 0))
})

test_that("hybrid curve follows KM to the cut-point, then is continuous", {
  fe <- fit_endpoint(example_ipd, "OS", example_arms)
  t <- make_grid(10, 1 / 52)
  s <- fit_curve(fe$indep$comp$weibull, "Standard care", t, fe$km$comp, TRUE, 1.5)
  before <- t <= 1.5
  expect_equal(s[before], km_surv(fe$km$comp, t[before]))
  i <- max(which(before))
  expect_lt(abs(s[i + 1] - s[i]), 0.01)              # no jump at the join
})

test_that("full pipeline runs and passes computational checks", {
  fits <- fit_all(example_ipd, example_arms)
  spec <- default_spec(lifetable = example_lt, fu_years = 4)
  out <- run_model(spec, fits)
  v <- out$validation
  expect_false(any(v$status[v$domain == "Computational"] == "fail"))
  expect_gt(out$res$inc_qaly, 0)
  # HR approach using the trial Cox HR
  cx <- cox_summary(fits$OS)
  spec$endpoints$OS$approach <- "hr"
  spec$endpoints$OS$hr <- list(hr = cx$hr, se = cx$se_log,
                               lower = exp(log(cx$hr) - 1.96 * cx$se_log),
                               upper = exp(log(cx$hr) + 1.96 * cx$se_log))
  out2 <- run_model(spec, fits)
  expect_equal(out2$curves$int$OS[1:50],
               (out2$curves$comp$OS^cx$hr)[1:50], tolerance = 1e-10)
})

test_that("DSA and PSA run and summarise", {
  fits <- fit_all(example_ipd, example_arms)
  spec <- default_spec(lifetable = example_lt, fu_years = 4)
  spec$endpoints$OS$approach <- "hr"
  spec$settings$cycle <- 1 / 12
  t <- make_grid(30, 1 / 12)
  pieces <- predict_pieces(spec, fits, t)
  d <- run_dsa(pieces, spec)
  expect_true(all(c("icer_low", "icer_high") %in% names(d)))
  runs <- run_psa(pieces, spec, n = 50, seed = 3)
  expect_equal(nrow(runs), 50)
  expect_true(all(runs$valid))
  expect_gt(stats::sd(runs$hr_os), 0)
})

test_that("flexsurv adapter matches survreg for shared models", {
  skip_if_not_installed("flexsurv")
  d <- example_ipd[example_ipd$endpoint == "OS" & example_ipd$arm == "Standard care", ]
  d$arm <- factor(d$arm)
  for (m in c("exp", "weibull", "lnorm", "llogis")) {
    a <- fit_one(d, m, FALSE, "flexsurv")
    b <- fit_one(d, m, FALSE, "survreg")
    expect_equal(a$loglik, b$loglik, tolerance = 1e-4, info = m)
    expect_equal(a$S(c(0.5, 2, 10)), b$S(c(0.5, 2, 10)), tolerance = 1e-4, info = m)
  }
  dd <- example_ipd[example_ipd$endpoint == "OS", ]
  dd$arm <- factor(dd$arm, levels = c("Standard care", "New drug"))
  a <- fit_one(dd, "weibull", TRUE, "flexsurv")
  b <- fit_one(dd, "weibull", TRUE, "survreg")
  expect_equal(a$S(c(1, 5), "New drug"), b$S(c(1, 5), "New drug"), tolerance = 1e-4)
  for (m in c("gompertz", "gamma", "gengamma", "spline1", "spline2")) {
    f <- fit_one(d, m, FALSE, "flexsurv")
    expect_true(f$ok, info = m)
    expect_true(all(diff(f$S(seq(0, 20, 0.25))) <= 1e-10), info = m)
  }
})
