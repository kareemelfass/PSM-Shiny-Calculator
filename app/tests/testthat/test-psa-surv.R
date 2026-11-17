# V2: survival-parameter uncertainty in the PSA

psa_setup <- function(engine = "survreg", cycle = 1 / 12) {
  fits <- list(OS = fit_endpoint(example_ipd, "OS", example_arms, engine),
               PFS = fit_endpoint(example_ipd, "PFS", example_arms, engine),
               ToT = fit_endpoint(example_ipd, "ToT", example_arms, engine))
  spec <- default_spec(example_arms, lifetable = example_lt, fu_years = 4)
  spec$settings$cycle <- cycle
  t <- make_grid(spec$settings$horizon, cycle)
  list(fits = fits, spec = spec, pieces = predict_pieces(spec, fits, t))
}

test_that("multivariate normal draws have the right mean and covariance", {
  set.seed(1)
  mu <- c(1, -2, 0.5)
  V <- matrix(c(0.04, 0.01, 0, 0.01, 0.09, -0.02, 0, -0.02, 0.01), 3)
  x <- mvn_draw(20000, mu, V)
  expect_equal(colMeans(x), mu, tolerance = 0.01)
  expect_equal(stats::cov(x), V, tolerance = 0.05)
  # Not positive definite (rank 1): eigen fallback still gives the right covariance
  V1 <- matrix(c(1, 1, 1, 1), 2)
  y <- mvn_draw(20000, c(0, 0), V1)
  expect_equal(stats::cov(y), V1, tolerance = 0.05)
})

test_that("survreg sampled curves are centred on the fitted curve", {
  d <- example_ipd[example_ipd$endpoint == "OS" & example_ipd$arm == "Standard care", ]
  d$arm <- factor(d$arm)
  f <- fit_one(d, "weibull", FALSE, "survreg")
  set.seed(2)
  draw <- f$draw_mvn(3000)
  s2 <- vapply(1:3000, \(i) draw(i, 2), numeric(1))
  expect_equal(stats::median(s2), f$S(2), tolerance = 0.01)
  expect_gt(stats::sd(s2), 0.005)
  expect_equal(draw(1, 0), 1)
})

test_that("joint model: both arms use the same draw (shared scale, PH holds per draw)", {
  d <- example_ipd[example_ipd$endpoint == "OS", ]
  d$arm <- factor(d$arm, levels = c("Standard care", "New drug"))
  f <- fit_one(d, "weibull", TRUE, "survreg")
  set.seed(3)
  draw <- f$draw_mvn(200)
  # Weibull PH: log(-log S_int) - log(-log S_comp) is constant over time within a draw.
  # Independent draws per arm would give different scales, so the gap would change with t.
  gap <- \(i, t) log(-log(draw(i, t, "New drug"))) - log(-log(draw(i, t, "Standard care")))
  diffs <- vapply(1:200, \(i) gap(i, 1) - gap(i, 6), numeric(1))
  expect_lt(max(abs(diffs)), 1e-8)
  # ...and the gap itself varies between draws (the treatment effect is uncertain)
  expect_gt(stats::sd(vapply(1:200, \(i) gap(i, 2), numeric(1))), 0.01)
})

test_that("flexsurv: parameter route reproduces every fitted curve", {
  skip_if_not_installed("flexsurv")
  d <- example_ipd[example_ipd$endpoint == "OS" & example_ipd$arm == "Standard care", ]
  d$arm <- factor(d$arm)
  for (m in available_models("flexsurv")$id) {
    f <- fit_one(d, m, FALSE, "flexsurv")
    expect_true(f$ok, info = m)
    set.seed(4)
    draw <- f$draw_mvn(200)                     # stops if the self-check fails
    s <- vapply(1:200, \(i) draw(i, 3), numeric(1))
    expect_true(all(s >= 0 & s <= 1), info = m)
    expect_equal(stats::median(s), f$S(3), tolerance = 0.03, info = m)
  }
  dd <- example_ipd[example_ipd$endpoint == "OS", ]
  dd$arm <- factor(dd$arm, levels = c("Standard care", "New drug"))
  for (m in c("weibull", "gengamma", "spline1")) {
    f <- fit_one(dd, m, TRUE, "flexsurv")
    expect_no_error(f$draw_mvn(50), message = m)
  }
})

test_that("flexsurv and survreg give the same Weibull uncertainty", {
  skip_if_not_installed("flexsurv")
  # Same model, different parameterisation (a linear map on the log scale),
  # so the sampled survival distribution must agree.
  d <- example_ipd[example_ipd$endpoint == "OS" & example_ipd$arm == "Standard care", ]
  d$arm <- factor(d$arm)
  a <- fit_one(d, "weibull", FALSE, "flexsurv")
  b <- fit_one(d, "weibull", FALSE, "survreg")
  set.seed(5)
  da <- a$draw_mvn(4000)
  db <- b$draw_mvn(4000)
  sa <- vapply(1:4000, \(i) da(i, 2), numeric(1))
  sb <- vapply(1:4000, \(i) db(i, 2), numeric(1))
  expect_equal(stats::sd(sa), stats::sd(sb), tolerance = 0.06)
  expect_equal(stats::median(sa), stats::median(sb), tolerance = 0.01)
})

test_that("bootstrap keeps patients whole and arm sizes fixed", {
  set.seed(6)
  b <- boot_ipd(example_ipd)
  n_orig <- table(example_ipd$arm, example_ipd$endpoint)
  n_boot <- table(b$arm, b$endpoint)
  expect_equal(n_boot, n_orig)
  # Every sampled patient brings all three endpoints
  per_id <- table(b$id) / table(b$id[b$endpoint == "OS"])[names(table(b$id))]
  expect_true(all(per_id == 3))
  # Without ids: resampled within arm x endpoint, sizes still fixed
  no_id <- example_ipd
  no_id$id <- NA
  b2 <- boot_ipd(no_id)
  expect_equal(table(b2$arm, b2$endpoint), n_orig)
})

test_that("fixed-curve parameters are identical whatever the survival method", {
  x <- psa_setup()
  r0 <- run_psa(x$pieces, x$spec, n = 30, seed = 7, surv_method = "none")
  r1 <- run_psa(x$pieces, x$spec, n = 30, seed = 7, surv_method = "mvn", fits = x$fits)
  expect_equal(r0[, c("hr_os", "hr_pfs", "u_pf", "u_pd")], r1[, c("hr_os", "hr_pfs", "u_pf", "u_pd")])
  # "none" reproduces the deterministic curves exactly in every run
  b <- attr(r0, "bands")
  expect_equal(unname(b$curves$comp_OS[1, ]), unname(b$curves$comp_OS[30, ]))
  expect_false(isTRUE(all.equal(b$curves$comp_OS, attr(r1, "bands")$curves$comp_OS)))
})

test_that("bootstrap PSA is reproducible and the same in serial and parallel", {
  x <- psa_setup()
  x$spec$tot$int$method <- "fit"        # also sample a ToT curve
  x$pieces <- predict_pieces(x$spec, x$fits, x$pieces$t)
  a <- run_psa(x$pieces, x$spec, n = 12, seed = 8, surv_method = "boot",
               fits = x$fits, ipd = example_ipd)
  b <- run_psa(x$pieces, x$spec, n = 12, seed = 8, surv_method = "boot",
               fits = x$fits, ipd = example_ipd)
  expect_equal(a$d_qaly, b$d_qaly)
  skip_if(parallel::detectCores() < 2, "needs 2 cores")
  p <- run_psa(x$pieces, x$spec, n = 12, seed = 8, surv_method = "boot",
               fits = x$fits, ipd = example_ipd, cores = 2, batch = 6,
               code_dir = file.path(app_dir, "R"))
  expect_true(all(p$valid), info = paste(unique(p$reason), collapse = "; "))
  expect_equal(p$d_qaly, a$d_qaly)
  expect_equal(p$cost_int, a$cost_int)
  # Same for parameter sampling (draws are made before the runs are split)
  m1 <- run_psa(x$pieces, x$spec, n = 12, seed = 8, surv_method = "mvn", fits = x$fits)
  m2 <- run_psa(x$pieces, x$spec, n = 12, seed = 8, surv_method = "mvn", fits = x$fits,
                cores = 2, batch = 6, code_dir = file.path(app_dir, "R"))
  expect_equal(m2$d_qaly, m1$d_qaly)
})

test_that("invalid runs are counted with a reason, never silently fixed", {
  x <- psa_setup()
  spec <- x$spec
  spec$settings$floor_on <- FALSE
  bad <- x$pieces
  bad$gp <- NULL
  bad$OS$comp[10:20] <- bad$OS$comp[10] + 0.01 * (1:11)   # OS rises
  draws <- list(hr_os = 1, hr_pfs = 1, u_pf = 0.75, u_pd = 0.6, run_seed = 1,
                costs = lapply(c(comp = "comp", int = "int"),
                               \(a) matrix(unlist(spec$econ[[a]][cost_fields]), 1,
                                           dimnames = list(NULL, cost_fields))))
  f <- make_run_fun(spec, \(i, seed) list(pieces = bad, fail = NULL), draws, 1:5)
  r <- f(1)$row
  expect_false(r$valid)
  expect_match(r$reason, "Survival increases")
  g <- make_run_fun(spec, \(i, seed) list(pieces = NULL, fail = "Refit failed: OS"), draws, 1:5)
  expect_match(g(1)$row$reason, "Refit failed")
  runs <- rbind(transform(f(1)$row, d_cost = 1, d_qaly = 1), transform(g(1)$row, d_cost = 1, d_qaly = 1))
  runs$valid[1] <- TRUE
  s <- psa_summary(runs, 30000)
  expect_equal(s$n_invalid, 1)
  expect_equal(names(s$reasons), "Refit failed: OS")
  # No valid runs at all: summary still works and shows why
  runs$valid <- FALSE
  s0 <- psa_summary(runs, 30000)
  expect_true(is.na(s0$icer))
  expect_equal(s0$n_invalid, 2)
})

test_that("summaries: bands, deterministic comparison, convergence", {
  x <- psa_setup()
  r <- run_psa(x$pieces, x$spec, n = 40, seed = 9, surv_method = "mvn", fits = x$fits)
  bd <- psa_bands(r, example_arms)
  expect_true(all(bd$lower <= bd$median + 1e-12 & bd$median <= bd$upper + 1e-12))
  det <- run_engine(assemble_curves(x$pieces, x$spec), x$spec)
  cmp <- psa_vs_deterministic(r, det)
  expect_equal(cmp$deterministic[cmp$measure == "d_qaly"], det$inc_qaly)
  expect_lt(abs(cmp$rel_diff[cmp$measure == "qaly_comp"]), 0.1)
  cv <- psa_convergence(r, 30000)
  expect_equal(nrow(cv), sum(r$valid))
  s <- psa_summary(r, 30000)
  expect_equal(tail(cv$icer, 1), s$icer)
})
