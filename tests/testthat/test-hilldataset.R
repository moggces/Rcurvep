context("model fit using a dataset")

data("zfishbeh")

test_that("use cc2", {

  outp <- run_fit(zfishbeh, modls = "cc2")
  expect_length(outp, 2)
  expect_true(all(outp$result$fit_set$win_modl == "cc2"))

})

test_that("output fit_set only", {
  outp <- run_fit(zfishbeh, keep_sets = "fit_set", modls = "hill")
  expect_length(outp$result, 1)
})

test_that("use only cnst model", {
  outp <- run_fit(zfishbeh, modls = "cnst")
  expect_true(all(outp$result$fit_set$win_modl == "cnst"))
})

test_that("set hill_pdir", {
  outp <- run_fit(zfishbeh, keep_sets = "fit_set", hill_pdir = -1, modls = "hill")
  tp_d <- outp$result$fit_set$hill_tp
  expect_true(all(na.omit(tp_d) < 0))
})

test_that("use n_samples ", {
  set.seed(300)
  outp <- run_fit(zfishbeh, n_samples = 2, modls = "hill")
  expect_true(all(outp$result$fit_set$sample_id %in% c(1, 2)))

})

test_that("warnings", {
  expect_warning(run_fit(zfishbeh, pdir = 1, modls = "cc2"))
})

test_that("errors", {
  expect_error(run_fit(data(zfishdev), modls = "cc2"))
  expect_error(run_fit(zfishbeh, modls = "xx"))
  expect_error(run_fit(zfishbeh, modls = "cc2", keep_sets = "resp_set"))
  expect_error(run_fit(zfishbeh, modls = "cc2", n_samples = 2))

})
