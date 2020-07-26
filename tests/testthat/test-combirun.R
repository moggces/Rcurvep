context("run rcurvep using a combination of parameters")

data("zfishbeh")
dat <- zfishbeh

test_that("no input parameters", {
  expect_error(combi_run_rcurvep(dat))
})

test_that("n_samples = NULL", {

  thresholds <- c(5, 10)
  ranges <- c(-100000, 100000)
  outp <- combi_run_rcurvep(dat, n_samples = NULL, TRSH = thresholds, RNGE = ranges, keep_sets = "act_set")
  expect_length(outp$result, 1)
  expect_equal(
    nrow(outp$result$act_set),
    nrow(dplyr::distinct(dat, chemical, endpoint))*length(thresholds)*length(ranges))
  expect_true(all(thresholds %in% outp$config$TRSH) && all(ranges %in% outp$config$RNGE) )

})

test_that("n_samples != NULL", {

  thresholds <- c(5, 10)
  n_samples <- 2
  outp <- combi_run_rcurvep(zfishbeh, n_samples = n_samples, TRSH = thresholds)
  expect_equal(
    nrow(outp$result$act_set),
    nrow(dplyr::distinct(dat, chemical, endpoint))*length(thresholds)*n_samples
  )
  expect_equal(dplyr::n_distinct(outp$result$act_set$TRSH), 2)
  expect_true(all(thresholds %in% outp$config$TRSH))
})



test_that("seed is not null", {

  thresholds <- c(5, 10)
  n_samples <- 2
  seed <- 300
  outp1 <- combi_run_rcurvep(zfishbeh, n_samples = n_samples, TRSH = thresholds, seed = seed)
  outp2 <- combi_run_rcurvep(zfishbeh, n_samples = n_samples, TRSH = thresholds, seed = seed)

  expect_true(identical(outp1$result$resp_set, outp2$result$resp_set))
  expect_true(outp1$config$seed == seed)
})

