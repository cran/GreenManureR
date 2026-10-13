toy <- data.frame(
  treatment  = rep(c("Control", "GM"), each = 4),
  block      = rep(1:4, 2),
  fert_N     = 0,
  gm_biomass = rep(c(0, 5), each = 4),
  gm_N       = rep(c(0, 2.5), each = 4),
  grain_yield = rep(c(3, 4), each = 4),
  grain_N    = 1,
  straw_yield = rep(c(4, 5), each = 4),
  straw_N    = 0.5)

test_that("N balance follows the documented formulas", {
  nb <- greenmanure_n_balance(toy, nfrv = FALSE)
  p <- nb$plot
  gm <- p[p$treatment == "GM", ][1, ]
  expect_equal(gm$N_gm_input, 5 * 1000 * 2.5 / 100)                 # 125
  expect_equal(gm$N_uptake, 4 * 1000 * 1 / 100 + 5 * 1000 * 0.5 / 100)  # 65
  expect_equal(gm$N_balance, 125 - 65)
  ctl_uptake <- 3 * 10 + 4 * 5                                       # 50
  expect_equal(gm$ANR, (65 - ctl_uptake) / 125 * 100)
  expect_equal(gm$AE_gmN, (4 - 3) * 1000 / 125)
  expect_true(all(is.na(p$ANR[p$treatment == "Control"])))
})

test_that("summary has one row per treatment and fertilizer level", {
  nb <- greenmanure_n_balance(gm_agronomy, by = "season")
  expect_equal(nrow(nb$summary), 2 * 4 * 3)
  expect_true(all(c("ANR_mean", "ANR_se") %in% names(nb$summary)))
})

test_that("NFRV is positive for green manure at zero fertilizer N", {
  nb <- greenmanure_n_balance(gm_agronomy, by = "season")
  s <- nb$summary
  k <- s[s$season == "Kharif" & s$fert_N == 0 & s$treatment != "Control", ]
  expect_true(all(k$NFRV_mean > 0))
})

test_that("missing control is an error", {
  expect_error(greenmanure_n_balance(toy, control = "None"), "not found")
})
