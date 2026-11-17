# Guyot reconstruction: round trip from simulated IPD

test_that("reconstruction recovers the KM curve and patient count", {
  set.seed(11)
  n <- 300
  ev <- rweibull(n, 1.3, 2.5)
  ce <- pmin(runif(n, 1, 6), 5)
  d <- data.frame(time = pmin(ev, ce), event = as.integer(ev <= ce))
  km <- survival::survfit(survival::Surv(time, event) ~ 1, data = d)
  digi <- data.frame(time = c(0, km$time), surv = c(1, km$surv))
  rt <- seq(0, 4.5, 0.5)
  risk <- data.frame(time = rt, n = summary(km, times = rt, extend = TRUE)$n.risk)

  rec <- guyot_reconstruct(digi, risk, tot_events = sum(d$event))
  expect_equal(nrow(rec), n)
  expect_equal(sum(rec$event), sum(d$event), tolerance = 3)

  km2 <- survival::survfit(survival::Surv(time, event) ~ 1, data = rec)
  tt <- c(0.5, 1, 2, 3, 4)
  expect_lt(max(abs(km_surv(km2, tt) - km_surv(km, tt))), 0.02)
})

test_that("example digitised file reconstructs close to the original", {
  digi <- read.csv(file.path(app_dir, "example", "example_digitised_os_standard_care.csv"))
  risk <- read.csv(file.path(app_dir, "example", "example_nrisk_os_standard_care.csv"))
  rec <- guyot_reconstruct(digi, risk)
  orig <- example_ipd[example_ipd$arm == "Standard care" & example_ipd$endpoint == "OS", ]
  k1 <- survival::survfit(survival::Surv(time, event) ~ 1, data = orig)
  k2 <- survival::survfit(survival::Surv(time, event) ~ 1, data = rec)
  # rec is in months (as digitised); orig is in years
  tt <- c(6, 12, 24, 36)
  expect_lt(max(abs(km_surv(k1, tt / 12) - km_surv(k2, tt))), 0.03)
})

test_that("bad input gives a clear error", {
  digi <- data.frame(time = c(0, 1, 2), surv = c(1, 0.8, 0.6))
  expect_error(guyot_reconstruct(digi, data.frame(time = c(1, 2), n = c(10, 5))), "start at time 0")
})
