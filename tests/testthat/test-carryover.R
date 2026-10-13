test_that("carryover_effect detects residual effects and retention", {
  co <- carryover_effect(gm_agronomy, strata = "fert_N",
                         season_levels = c("Kharif", "Rabi"))
  expect_s3_class(co, "gm_carryover")
  e <- co$effects
  expect_equal(nrow(e), 3 * 2 * 3)
  expect_true(all(e$diff[e$season == "Kharif" & e$fert_N == 0] > 0))
  expect_equal(e$retention_vs_first[e$season == "Kharif"],
               rep(100, sum(e$season == "Kharif")))
  # residual effect in the second crop is smaller than in the first
  expect_true(all(e$retention_vs_first[e$season == "Rabi"] < 100))
  expect_true(all(e$p_adj >= e$p - 1e-12))
  expect_false(is.null(co$pooled))
  expect_equal(nrow(co$cumulative), 3 * 3)
})

test_that("contrasts equal mean differences in a balanced design", {
  d <- gm_agronomy[gm_agronomy$fert_N == 0 & gm_agronomy$season == "Kharif", ]
  co <- carryover_effect(d, season_levels = "Kharif")
  m <- tapply(d$grain_yield, d$treatment, mean)
  e <- co$effects[co$effects$treatment == "Mucuna", ]
  expect_equal(e$diff, unname(m["Mucuna"] - m["Control"]), tolerance = 1e-10)
  expect_null(co$pooled)
})

test_that("carryover_effect validates input", {
  expect_error(carryover_effect(gm_agronomy, control = "None"), "not found")
  expect_error(carryover_effect(gm_agronomy, response = "zzz"), "not found")
})

test_that("joint analysis links the four analyses", {
  fit <- decomp_fit(gm_litterbag, group = "treatment")
  nr  <- nutrient_release(gm_litterbag, group = "treatment", biomass = "biomass_tha")
  nb  <- greenmanure_n_balance(gm_agronomy, by = "season")
  co  <- carryover_effect(gm_agronomy, strata = "fert_N")
  expect_warning(ja <- gm_joint_analysis(fit, nr, nb, co, fert_N_level = 0),
                 "descriptive")
  expect_s3_class(ja, "gm_joint")
  expect_equal(nrow(ja$data), 3)
  expect_true(all(c("t50", "N_release_d60", "ANR", "yield_gain_pct_Kharif") %in% names(ja$data)))
  expect_error(gm_joint_analysis(fit, nr, nb, co, nutrient_name = "Zn"), "not found")
})
