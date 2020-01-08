context("model fit using a dataset")

data("zfishbeh")

test_that("use defaults", {

  outp <- run_fit(zfishbeh)
  expect_length(outp, 2)

})

test_that("output fit_set only", {
  outp <- run_fit(zfishbeh, keep_sets = "fit_set")
  expect_length(outp$result, 1)
})

test_that("use only cnst model", {
  outp <- run_fit(zfishbeh, modls = "cnst")
  expect_true(all(outp$result$fit_set$win_modl == "cnst"))
})

test_that("set hill_pdir", {
  outp <- run_fit(zfishbeh, keep_sets = "fit_set", hill_pdir = -1)
  tp_d <- outp$result$fit_set$hill_tp
  expect_true(all(na.omit(tp_d) < 0))
})

test_that("use n_samples ", {
  set.seed(300)
  outp <- run_fit(zfishbeh, n_samples = 2)
  expect_true(all(outp$result$fit_set$sample_id %in% c(1, 2)))

})

test_that("warnings", {
  expect_warning(run_fit(zfishbeh, pdir = 1))
})

test_that("errors", {
  expect_error(run_fit(data(zfishdev)))
  expect_error(run_fit(zfishbeh, modls = "xx"))
  expect_error(run_fit(zfishbeh, keep_sets = "resp_set"))

})
