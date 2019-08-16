context("run rcurvep using a combination of parameters")

data("zfishbeh")
dat <- create_dataset(zfishbeh)


test_that("no input parameters", {
  expect_error(combi_run_rcurvep(dat))
})

test_that("n_samples = NULL", {

  thresholds <- c(5, 10)
  ranges <- c(-100000, 100000)
  outp <- combi_run_rcurvep(dat, n_samples = NULL,
                            keep = c("act_set", "resp_set"), TRSH = thresholds, RNGE = ranges)
  expect_length(outp$result, 2)
  expect_equal(
    nrow(outp$result$act_set),
    nrow(dplyr::distinct(dat, chemical, endpoint))*length(thresholds)*length(ranges))

})

test_that("n_samples != NULL", {

  thresholds <- c(5, 10)
  n_samples <- 2
  outp <- combi_run_rcurvep(zfishbeh, n_samples = n_samples, keep = c("act_set"), TRSH = thresholds)
  expect_equal(
    nrow(outp$result$act_set),
    nrow(dplyr::distinct(dat, chemical, endpoint))*length(thresholds)*n_samples
  )

})


test_that("on a single rcurvep object", {

  outp <- run_rcurvep(dat)
  obj <- merge_rcurvep_output(outp, keep_data = "act_set")
  expect_equal(
    nrow(obj$result$act_set),
    nrow(dplyr::distinct(dat, chemical, endpoint))
  )

})
