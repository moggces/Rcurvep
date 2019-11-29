context("basic curvep running")

# make some datasets
data("zfishbeh")
dats <- create_dataset(zfishbeh) %>%
  split(.$endpoint)

test_that("run one resp per conc per endpoint", {

  outp <- run_rcurvep(dats[[1]], config = curvep_defaults())
  expect_length(outp, 2)

})

test_that("run one resp per conc many endpoints", {
  outp <- run_rcurvep(dplyr::bind_rows(dats), config = curvep_defaults())
  expect_length(outp, 2)
})

test_that("with parameter change", {
  outp <- run_rcurvep(dats[[1]], TRSH = 30)
  expect_equal(outp$config$TRSH, 30)
})


test_that("with mask", {
  outp <- run_rcurvep(dats[[1]], mask = c(1, 2))
  expect_equal(sum(tail(outp$result$resp_set$mask, n = 5)), 2)
})

test_that("only one set", {
  outp <- run_rcurvep(dats[[1]], keep_sets = "act_set")
  expect_true(length(outp$result) == 1 && ('act_set' %in% names(outp$result)))
})

test_that("dataset has mask column", {
  inp <- dats[[1]] %>%
    dplyr::mutate(mask = rep(c(0, 0, 0, 0, 1),2))

  outp <- run_rcurvep(inp)
  expect_equal(sum(tail(outp$result$resp_set$mask, n = 5)), 1)
})





