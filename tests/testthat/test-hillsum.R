context("summary of model fit of a dataset")

data("zfishbeh")
fitd1 <- run_fit(zfishbeh, modls = "cc2")
fitd2 <- run_fit(zfishbeh, n_samples = 3, modls = "hill")
fitd3 <- run_fit(zfishbeh, modls = "cnst")


test_that("use defaults", {

  expect_length(summarize_fit_output(fitd1), 3)
  expect_length(summarize_fit_output(fitd2), 3)
  expect_length(summarize_fit_output(fitd3), 3)

})
