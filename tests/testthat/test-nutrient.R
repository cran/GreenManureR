test_that("nutrient_release computes remaining nutrient and kg/ha", {
  nr <- nutrient_release(gm_litterbag, nutrients = c("N", "K"),
                         group = "treatment", biomass = "biomass_tha",
                         at = c(30, 60))
  expect_s3_class(nr, "gm_nutrient")
  d0 <- nr$data[nr$data$time == 0 & nr$data$nutrient == "N", ]
  expect_equal(mean(d0$remaining), 100, tolerance = 1e-8)
  expect_true(all(c("input_kg_ha", "released_kg_ha") %in% names(nr$data)))
  # K is leached faster than N
  t50 <- nr$t50
  expect_true(all(t50$t50[t50$nutrient == "K"] < t50$t50[t50$nutrient == "N"]))
  expect_equal(nrow(nr$release_at), 3 * 2 * 2)
})

test_that("immobilization gives remaining nutrient above 100%", {
  d <- data.frame(time = rep(c(0, 7, 14, 28, 42, 56), each = 2),
                  mass = rep(c(100, 90, 80, 70, 60, 50), each = 2))
  d$P <- 0.2 * (1 + 0.5 * (d$time / 56)) * 100 / d$mass
  nr <- nutrient_release(d, nutrients = "P", models = "single")
  expect_gt(max(nr$data$remaining), 100)
})

test_that("nutrient_release validates input", {
  expect_error(nutrient_release(gm_litterbag, nutrients = "Zn"), "not found")
})
