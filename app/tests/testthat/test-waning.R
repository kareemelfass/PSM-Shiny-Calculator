# V3 waning: reproduces Doc 12's worked example and the new options

# Doc 12 example: exponential OS, comparator hazard 0.18/yr, intervention 0.09/yr (HR 0.5),
# effect ceases at year 3.
doc_t <- make_grid(15, 1 / 365.25)
doc_c <- exp(-0.18 * doc_t)
doc_i <- exp(-0.09 * doc_t)
at <- \(s, y) stats::approx(doc_t, s, xout = y)$y
pts <- \(x) round(100 * x, 1)

test_that("Doc 12 numbers: survival at years 3 and 10 under each assumption", {
  expect_equal(pts(at(doc_i, 3)), 76.3)
  expect_equal(pts(at(doc_c, 3)), 58.3)
  expect_equal(pts(at(doc_c, 10)), 16.5)                                     # comparator
  expect_equal(pts(at(doc_i, 10)), 40.7)                                     # sustained effect
  eq <- apply_waning(doc_i, doc_c, doc_t, "equal", start = 3)
  expect_equal(pts(at(eq, 10)), 21.7)                                        # equal future hazards
  sd <- apply_waning(doc_i, doc_c, doc_t, "stop_drop", start = 3)
  expect_equal(pts(at(sd, 10)), 16.5)                                        # stop and drop
})

test_that("Doc 12 numbers: the gained advantage persists but narrows (18.1, 8.8, 5.1 points)", {
  eq <- apply_waning(doc_i, doc_c, doc_t, "equal", start = 3)
  expect_equal(pts(at(eq, 3) - at(doc_c, 3)), 18.1)
  expect_equal(pts(at(eq, 7) - at(doc_c, 7)), 8.8)
  expect_equal(pts(at(eq, 10) - at(doc_c, 10)), 5.1)
  # Equal hazards do not mean equal survival
  expect_true(all(eq[doc_t > 3] > doc_c[doc_t > 3]))
})

test_that("defaults (target 1, no frailty) give exactly the V2 result", {
  v2 <- \(s_i, s_c, start, end) {
    w <- waning_weight(cycle_mid(doc_t), start, end)
    from_hazard((1 - w) * cycle_hazard(s_i) + w * cycle_hazard(s_c))
  }
  expect_equal(apply_waning(doc_i, doc_c, doc_t, "gradual", 3, 5), v2(doc_i, doc_c, 3, 5), tolerance = 1e-12)
  expect_equal(apply_waning(doc_i, doc_c, doc_t, "equal", 3), v2(doc_i, doc_c, 3, 3), tolerance = 1e-12)
})

test_that("partial waning reaches the target HR, not 1", {
  w <- apply_waning(doc_i, doc_c, doc_t, "gradual", 3, 5, target_hr = 0.8)
  hr <- cycle_hazard(w) / cycle_hazard(doc_c)
  tm <- cycle_mid(doc_t)
  expect_equal(hr[tm < 3][1], 0.5, tolerance = 1e-8)
  expect_equal(hr[which.min(abs(tm - 4))], 0.65, tolerance = 1e-3)          # half way
  expect_equal(utils::tail(hr, 1), 0.8, tolerance = 1e-8)
})

test_that("frailty: individual-level waning pushes the marginal HR above 1 (Jennings 2024)", {
  theta <- 0.5
  w <- apply_waning(doc_i, doc_c, doc_t, "equal", start = 3, frailty = theta)
  hr <- cycle_hazard(w) / cycle_hazard(doc_c)
  tm <- cycle_mid(doc_t)
  # Analytic: marginal HR just after t0 = (1 + theta Hc_comp(t0)) / (1 + theta Hc_int(t0)),
  # with conditional cumulative hazards Hc = (exp(theta Hm) - 1) / theta
  hc <- \(Hm) expm1(theta * Hm) / theta
  expected <- (1 + theta * hc(0.18 * 3)) / (1 + theta * hc(0.09 * 3))
  expect_equal(hr[which(tm > 3)[1]], expected, tolerance = 1e-4)
  expect_gt(expected, 1.1)
  # So the benefit is smaller than under equal marginal hazards
  eq_marg <- apply_waning(doc_i, doc_c, doc_t, "equal", start = 3)
  expect_lt(at(w, 10), at(eq_marg, 10))
  # Before waning nothing changes; a tiny frailty matches the marginal version
  expect_equal(w[doc_t < 3], doc_i[doc_t < 3], tolerance = 1e-10)
  tiny <- apply_waning(doc_i, doc_c, doc_t, "equal", start = 3, frailty = 1e-8)
  expect_equal(tiny, eq_marg, tolerance = 1e-6)
})

test_that("timing can be anchored to the end of treatment (TA650)", {
  spec <- default_spec()
  spec$tot$int$cap <- 2
  spec$waning <- modifyList(spec$waning, list(type = "equal", anchor = "tx_stop", start = 3))
  wt <- waning_times(spec)
  expect_equal(wt$start, 5)                       # 3 years after stopping = 5 years after start
  expect_true(wt$anchored_to_cap)
  spec$tot$int$cap <- NA
  wt2 <- waning_times(spec)
  expect_equal(wt2$start, 3)
  expect_true(wt2$anchor_missing_cap)
  # The assembled curves and their diagnostics carry the absolute times
  spec$tot$int$cap <- 2
  p <- list(t = doc_t, OS = list(comp = doc_c, int = doc_i), PFS = list(comp = doc_c^2, int = doc_i^2),
            ToT = list(comp = NULL, int = NULL), gp = NULL)
  cv <- assemble_curves(p, spec)
  expect_equal(cv$waning$start, 5)
  expect_equal(waning_diagnostics(cv, spec)$start, 5)
  expect_equal(cv$int$OS[doc_t < 5], doc_i[doc_t < 5])          # full effect until year 5
})

test_that("diagnostics: stop-and-drop jump size and post-waning HR", {
  spec <- default_spec()
  t <- doc_t
  p <- list(t = t, OS = list(comp = doc_c, int = doc_i), PFS = list(comp = doc_c^2, int = doc_i^2),
            ToT = list(comp = NULL, int = NULL), gp = NULL)
  spec$waning <- modifyList(spec$waning, list(type = "stop_drop", start = 3))
  cv <- assemble_curves(p, spec)
  dg <- waning_diagnostics(cv, spec)
  expect_equal(pts(dg$OS$drop_pts), 18.1, tolerance = 0.11)
  spec$waning$type <- "equal"
  dg2 <- waning_diagnostics(assemble_curves(p, spec), spec)
  expect_equal(dg2$OS$hr_after, 1, tolerance = 1e-6)
  expect_equal(dg2$OS$share_above_target, 0)
  # Crossing hazards in the tail: flagged, but weighted by patients at risk
  pc <- p
  pc$OS$int <- exp(-(0.05 * doc_t + 0.02 * doc_t^2))      # intervention hazard rises above comparator late
  dg3 <- waning_diagnostics(assemble_curves(pc, spec), spec)
  expect_gt(dg3$OS$share_above_target, 0)
  expect_lt(dg3$OS$share_above_weighted, dg3$OS$share_above_target)
  expect_equal(dg3$OS$first_above, (0.18 - 0.05) / 0.04, tolerance = 0.01)   # hazards cross at 3.25 years
})

test_that("explorer, DSA and PSA include waning timing", {
  fits <- fit_all(example_ipd, example_arms, "survreg")
  spec <- default_spec(example_arms, lifetable = example_lt, fu_years = 4)
  spec$endpoints$OS$approach <- "hr"
  spec$settings$cycle <- 1 / 12
  pieces <- predict_pieces(spec, fits, make_grid(30, 1 / 12))

  g <- waning_grid(pieces, spec, starts = c(3, 5))
  expect_equal(nrow(g), 1 + 3 * 2)
  base_none <- run_engine(assemble_curves(pieces, spec), spec)
  expect_equal(g$icer[g$type == "none"], base_none$icer)
  # More waning, less benefit: equal-from-3 below equal-from-5; stop-and-drop lowest
  q <- \(ty, st) g$inc_qaly[g$type == ty & g$start_input %in% st]
  expect_lt(q("equal", 3), q("equal", 5))
  expect_lt(q("stop_drop", 5), q("equal", 5))

  expect_false(any(grepl("Waning", vapply(dsa_params(spec), `[[`, "", "name"))))
  spec$waning <- modifyList(spec$waning, list(type = "gradual", start = 5, end = 8))
  nm <- vapply(dsa_params(spec), `[[`, "", "name")
  expect_true(all(c("Waning start (years, end moves with it)", "Waning duration (years)") %in% nm))
  d <- run_dsa(pieces, spec)
  w <- d[d$parameter == "Waning start (years, end moves with it)", ]
  expect_lt(w$nmb_low, w$nmb_high)              # earlier waning, less benefit

  r0 <- run_psa(pieces, spec, n = 20, seed = 3)
  r1 <- run_psa(pieces, spec, n = 20, seed = 3, wane_start_range = c(3, 7))
  expect_true(all(r1$wane_start >= 3 & r1$wane_start <= 7))
  expect_true(all(is.na(r0$wane_start)))
  expect_equal(r0$hr_os, r1$hr_os)              # other draws unchanged
  expect_false(isTRUE(all.equal(r0$d_qaly, r1$d_qaly)))
})
