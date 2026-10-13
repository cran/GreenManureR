make_single <- function(k = 0.03, times = c(0, 7, 14, 28, 42, 56, 84, 112)) {
  data.frame(time = times, mass = 100 * exp(-k * times))
}

test_that("decomp_fit recovers a known single-pool rate", {
  fit <- decomp_fit(make_single(0.03), models = "single")
  expect_s3_class(fit, "gm_decomp")
  expect_equal(fit$table$k, 0.03, tolerance = 1e-4)
  expect_true(fit$table$selected)
})

test_that("decomp_fit selects the double-pool model for two-pool data", {
  t <- c(0, 7, 14, 21, 28, 42, 56, 84, 112, 140)
  d <- data.frame(time = t,
                  mass = 100 * (0.5 * exp(-0.08 * t) + 0.5 * exp(-0.004 * t)))
  fit <- decomp_fit(d)
  best <- fit$table[fit$table$selected, ]
  expect_equal(best$model, "double")
  expect_equal(best$k, 0.08, tolerance = 0.02)
  expect_gt(best$k, best$k_slow)
})

test_that("decomp_fit handles groups and mass on a fraction scale", {
  d <- rbind(cbind(make_single(0.02), sp = "a"), cbind(make_single(0.05), sp = "b"))
  d$mass <- d$mass / 100
  fit <- decomp_fit(d, group = "sp", models = "single", mass_scale = "fraction")
  expect_equal(nrow(fit$table), 2)
  expect_equal(fit$table$k[fit$table$sp == "b"], 0.05, tolerance = 1e-4)
})

test_that("decomp_fit validates input", {
  expect_error(decomp_fit(make_single(), time = "nope"), "not found")
  expect_error(decomp_fit(make_single()[1:2, ]), "No model could be fitted")
})

test_that("decomp_t50 matches log(2)/k for the single-pool model", {
  fit <- decomp_fit(make_single(0.03), models = "single")
  t50 <- decomp_t50(fit, ci = FALSE)
  expect_equal(t50$t50, log(2) / 0.03, tolerance = 1e-4)
  t95 <- decomp_t50(fit, p = 0.95, ci = FALSE)
  expect_true("t95" %in% names(t95))
  expect_equal(t95$t95, -log(0.05) / 0.03, tolerance = 1e-4)
})

test_that("decomp_t50 intervals bracket the estimate and are reproducible", {
  set.seed(10)
  d <- make_single(0.03); d$mass <- d$mass * (1 + rnorm(nrow(d), 0, 0.03))
  d$mass[1] <- 100
  fit <- decomp_fit(d, models = "single")
  a <- decomp_t50(fit)
  b <- decomp_t50(fit)
  expect_equal(a, b)
  expect_lt(a$lwr, a$t50); expect_gt(a$upr, a$t50)
})

test_that("decomp_t50 does not change the global RNG state", {
  fit <- decomp_fit(make_single(), models = "single")
  set.seed(1); x <- runif(1)
  set.seed(1); invisible(decomp_t50(fit)); y <- runif(1)
  expect_equal(x, y)
})

test_that("decomp_t50 validates its arguments", {
  fit <- decomp_fit(make_single(), models = "single")
  expect_error(decomp_t50(fit, p = 1.2), "strictly between")
  expect_error(decomp_t50(data.frame()), "decomp_fit")
})

test_that("asymptotic model returns NA when the target is never reached", {
  t <- c(0, 7, 14, 28, 42, 56, 84, 112)
  d <- data.frame(time = t, mass = 60 + 40 * exp(-0.05 * t))
  fit <- decomp_fit(d, models = "asymptotic")
  expect_true(is.na(decomp_t50(fit, p = 0.5, ci = FALSE)$t50))
})

test_that("predict returns model predictions", {
  fit <- decomp_fit(make_single(0.03), models = "single")
  pr <- predict(fit, times = c(0, 30))
  expect_equal(pr$pred, 100 * exp(-0.03 * c(0, 30)), tolerance = 1e-4)
})
