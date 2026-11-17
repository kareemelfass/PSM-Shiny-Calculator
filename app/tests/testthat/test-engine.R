# Known-answer tests for the economic engine and curve assembly

zero_econ <- function(spec) {
  for (a in c("comp", "int")) for (f in names(spec$econ[[a]])) {
    if (!f %in% c("interval_days", "rdi", "prog_share")) spec$econ[[a]][[f]] <- 0
  }
  spec
}

test_that("life-years match the exponential closed form", {
  spec <- default_spec()
  spec$settings$horizon <- 20
  spec$settings$cycle <- 1 / 365.25
  spec$settings$dr_out <- 0
  spec$settings$dr_cost <- 0
  t <- make_grid(20, spec$settings$cycle)
  cv <- assemble_curves(exp_pieces(spec, t), spec)
  r <- run_engine(cv, spec)
  expect_equal(r$comp$summary[["ly_undisc"]], (1 - exp(-0.3 * 20)) / 0.3, tolerance = 1e-5)
  expect_equal(r$int$summary[["ly_undisc"]], (1 - exp(-0.2 * 20)) / 0.2, tolerance = 1e-5)
  # PF life-years = area under PFS
  expect_equal(r$comp$summary[["ly_pf"]], (1 - exp(-0.8 * 20)) / 0.8, tolerance = 1e-5)
})

test_that("occupancy sums to one and stays within 0-1", {
  spec <- default_spec()
  t <- make_grid(30, 7 / 365.25)
  cv <- assemble_curves(exp_pieces(spec, t), spec)
  expect_true(occupancy_ok(cv))
})

test_that("HR = 1 reproduces the reference curve; S_int = S_comp^HR otherwise", {
  spec <- default_spec()
  spec$endpoints$OS$approach <- "hr"
  t <- make_grid(30, 7 / 365.25)
  p <- exp_pieces(spec, t)
  cv1 <- assemble_curves(p, spec, hr = list(OS = 1, PFS = 1))
  expect_equal(cv1$int$OS, cv1$comp$OS)
  cv2 <- assemble_curves(p, spec, hr = list(OS = 0.7, PFS = 1))
  expect_equal(cv2$int$OS, cv2$comp$OS^0.7)
  # Doc 9 worked example: S_B(24 months) = 0.40, HR 0.824 -> 0.470
  expect_equal(0.40^0.824, 0.470, tolerance = 1e-3)
})

test_that("HR direction is inverted correctly (Doc 9 B4a)", {
  h <- hr_to_log(1 / 0.7, 1 / 1.02, 1 / 0.48, "comp_vs_int")
  expect_equal(h$hr, 0.7)
  expect_equal(c(h$lower, h$upper), c(0.48, 1.02))
})

test_that("gradual waning keeps the survival already gained (Doc 12)", {
  t <- make_grid(15, 1 / 52)
  s_c <- exp(-0.18 * t)
  s_i <- exp(-0.09 * t)
  w <- apply_waning(s_i, s_c, t, "equal", start = 3)
  # After year 3 hazards are equal: the ratio of survival stays constant
  after <- t >= 3.1
  ratio <- w[after] / s_c[after]
  expect_lt(max(ratio) - min(ratio), 1e-3)
  # Gap at year 3 = exp(-0.27) - exp(-0.54) = 18.1 points (Doc 12)
  expect_equal(exp(-0.27) - exp(-0.54), 0.181, tolerance = 1e-3)
  # Gap at year 10 under equal hazards: 5.1 points (Doc 12)
  gap10 <- stats::approx(t, w - s_c, xout = 10)$y
  expect_equal(gap10, 0.051, tolerance = 2e-3)
  # Stop-and-drop equals the comparator from the start time
  sd <- apply_waning(s_i, s_c, t, "stop_drop", start = 3)
  expect_equal(sd[t >= 3], s_c[t >= 3])
})

test_that("gradual waning moves the HR linearly to 1 (Doc 9 C4a)", {
  t <- make_grid(15, 1 / 52)
  s_c <- exp(-0.2 * t)
  w <- apply_waning(s_c^0.72, s_c, t, "gradual", start = 3, end = 10)
  hr_t <- cycle_hazard(w) / cycle_hazard(s_c)
  tm <- cycle_mid(t)
  expect_equal(hr_t[tm < 3][1], 0.72, tolerance = 1e-8)
  expect_equal(hr_t[which.min(abs(tm - 6.5))], 0.86, tolerance = 2e-3)
  expect_equal(tail(hr_t, 1), 1, tolerance = 1e-8)
})

test_that("general-population floor binds and is reported", {
  spec <- default_spec(lifetable = example_lt)
  spec$settings$age <- 70
  t <- make_grid(40, 1 / 12)
  p <- exp_pieces(spec, t, lam_os = 0.01, lam_os_int = 0.01, lam_pfs = 0.02, lam_pfs_int = 0.02)
  p$gp <- gen_pop_surv(t, example_lt, 70, 0.5)
  cv <- assemble_curves(p, spec)
  expect_false(is.na(cv$flags$comp$floor_from))
  expect_true(all(cycle_hazard(cv$comp$OS) >= cycle_hazard(p$gp) - 1e-12))
  expect_true(all(cv$comp$OS <= p$gp + 1e-12))
})

test_that("PFS above OS is detected, constrained only when asked, and reported", {
  spec <- default_spec()
  t <- make_grid(10, 1 / 12)
  p <- exp_pieces(spec, t, lam_os = 0.5, lam_pfs = 0.3)   # PFS above OS: invalid
  cv <- assemble_curves(p, spec)
  expect_gt(cv$flags$comp$pfs_gt_os_cycles, 0)
  expect_true(occupancy_ok(cv))
  spec$settings$pfs_le_os <- FALSE
  cv2 <- assemble_curves(p, spec)
  expect_false(occupancy_ok(cv2))
})

test_that("drug cost follows time on treatment, with the cap applied", {
  spec <- zero_econ(default_spec())
  spec$settings$dr_cost <- 0
  spec$settings$horizon <- 10
  spec$settings$cycle <- 1 / 365.25
  spec$econ$int$drug_cost <- 1000
  spec$econ$int$interval_days <- 365.25     # 1000 per year on treatment
  spec$econ$int$rdi <- 1
  spec$tot$int <- list(method = "pfs", model = "weibull", cap = 2, beyond_prog = FALSE)
  t <- make_grid(10, spec$settings$cycle)
  cv <- assemble_curves(exp_pieces(spec, t), spec)
  r <- run_engine(cv, spec)
  # ToT = PFS (rate 0.5) capped at 2 years: area = (1 - exp(-1)) / 0.5
  expect_equal(r$int$summary[["cost_drug"]], 1000 * (1 - exp(-1)) / 0.5, tolerance = 1e-3)
})

test_that("ICER labels and ratio of means (Doc 13 Case 17)", {
  expect_match(icer_label(-10, 1), "Dominant")
  expect_match(icer_label(10, -1), "Dominated")
  runs <- data.frame(d_cost = c(5000, 10000, 15000), d_qaly = c(0.1, 0.2, 1.5),
                     valid = TRUE, constraint_active = FALSE)
  s <- psa_summary(runs, 30000)
  expect_equal(s$icer, 10000 / 0.6)       # ratio of means = 16,667
  expect_false(isTRUE(all.equal(s$icer, mean(runs$d_cost / runs$d_qaly))))
})

test_that("discounting lowers outcomes and costs", {
  spec <- default_spec()
  t <- make_grid(30, 7 / 365.25)
  cv <- assemble_curves(exp_pieces(spec, t), spec)
  r <- run_engine(cv, spec)
  expect_lt(r$int$summary[["ly_disc"]], r$int$summary[["ly_undisc"]])
  expect_lt(r$int$summary[["qaly"]], r$int$summary[["qaly_undisc"]])
})
