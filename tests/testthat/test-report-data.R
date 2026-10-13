test_that("simulate_greenmanure is reproducible and leaves the RNG alone", {
  a <- simulate_greenmanure(seed = 3); b <- simulate_greenmanure(seed = 3)
  expect_equal(a, b)
  set.seed(1); x <- runif(1)
  set.seed(1); invisible(simulate_greenmanure(seed = 9)); y <- runif(1)
  expect_equal(x, y)
})

test_that("shipped data sets match the simulator", {
  sim <- simulate_greenmanure(seed = 42)
  expect_equal(gm_litterbag, sim$litterbag)
  expect_equal(gm_agronomy, sim$agronomy)
})

test_that("print methods run and return the object invisibly", {
  fit <- decomp_fit(gm_litterbag, group = "treatment")
  expect_output(res <- print(fit), "gm_decomp")
  expect_identical(res, fit)
  nr <- nutrient_release(gm_litterbag, group = "treatment")
  expect_output(print(nr), "gm_nutrient")
})

test_that("greenmanure_report requires outdir and writes only there", {
  fit <- decomp_fit(gm_litterbag, group = "treatment")
  expect_error(greenmanure_report(decomp = fit), "outdir")
  skip_if_not_installed("ggplot2")
  out <- file.path(tempfile("gmrep"))
  expect_message(res <- greenmanure_report(decomp = fit, outdir = out, docx = FALSE),
                 "written")
  expect_true(file.exists(file.path(out, "T1_decomposition_models.csv")))
  expect_true(file.exists(file.path(out, "F1_decomposition_curves.png")))
  unlink(out, recursive = TRUE)
})

test_that("plot methods return ggplot objects", {
  skip_if_not_installed("ggplot2")
  fit <- decomp_fit(gm_litterbag, group = "treatment")
  expect_s3_class(plot(fit), "ggplot")
  nb <- greenmanure_n_balance(gm_agronomy, by = "season")
  expect_s3_class(plot(nb), "ggplot")
  co <- carryover_effect(gm_agronomy, strata = "fert_N")
  expect_s3_class(plot(co), "ggplot")
})
