context("basic curvep running")

# make some datasets
data("zfishbeh")
dats <- zfishbeh %>%
  dplyr::group_by(endpoint, chemical, conc) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
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
  expect_equal(sum(tail(outp$result[1,]$input[[1]]$mask, n = 2)), 2)
})

test_that("reparam by using input column from output", {
  outp <- run_rcurvep(dats[[1]])
  inp <- outp$result %>%
    dplyr::select(endpoint, chemical, input) %>%
    tidyr::unnest()
  outp2 <- run_rcurvep(inp, TRSH = 30)
  expect_equal(outp2$config$TRSH, 30)
})

test_that("reparam by using input column from output and new mask", {
  outp <- run_rcurvep(dats[[1]])
  inp <- outp$result %>%
    dplyr::select(endpoint, chemical, input) %>%
    tidyr::unnest()
  outp2 <- run_rcurvep(inp, mask = c(1, 2))
  expect_equal(sum(tail(outp2$result[1,]$input[[1]]$mask, n = 2)), 2)
})


# test_that("percent data with n_sample = NULL", {
#   data("zfishdev")
#   x <- zfishdev %>%
#     split(.$endpoint)
#   outd <- run_curvep_batch(x[[1]],
#                          directionality = 1,
#                          n_sample = NULL,
#                          threshold = 15,
#                          other_paras = list(CARR = 20, TrustHi = TRUE))
#   expect_true(sum(colnames(outd) %in% c('input', 'output', 'activity')) == 3, info = "a failed run" )
# })
#
#
# test_that("resp data with n_sample = NULL", {
#   data("zfishbeh")
#   x <- zfishbeh %>%
#     split(.$endpoint)
#   outd <- run_curvep_batch(x[[1]],
#                          directionality = 1,
#                          n_sample = NULL,
#                          threshold = 15,
#                          other_paras = list(CARR = 20, TrustHi = TRUE))
#   expect_true(sum(colnames(outd) %in% c('input', 'output', 'activity')) == 3, info = "a failed run" )
# })
#
#
# test_that("resp data original with multiple thresholds", {
#   data("zfishbeh")
#   x <- zfishbeh %>%
#     split(.$endpoint)
#   outd <- run_curvep_batch(x[[1]],
#                          directionality = 0,
#                          n_sample = NULL,
#                          threshold = list("1" = c(15,20), "-1" = c(30)),
#                          other_paras = list(CARR = 20, TrustHi = TRUE))
#
#   outd <- outd %>% dplyr::select(threshold, direction) %>% dplyr::distinct() %>%
#     dplyr::arrange(direction, threshold) %>% dplyr::pull(threshold)
#   expect_equal(outd, c(30, 15, 20))
# })
#
#
# test_that("percent boot dataset job run", {
#   data("zfishdev")
#   x <- zfishdev %>%
#     split(.$endpoint)
#   outd <- run_curvep_batch(x[[1]],
#                          directionality = 1,
#                          n_sample = 2,
#                          threshold = seq(5, 10, by = 5),
#                          other_paras = list(CARR = 20, TrustHi = TRUE))
#
#   outd <- outd %>% dplyr::select(threshold, direction, repeat_id) %>% dplyr::distinct()
#
#   expect_true(nrow(outd) == 4, info = "a failed run" )
# })
#
#
# test_that("resp boot dataset job run", {
#   data("zfishbeh")
#   x <- zfishbeh %>%
#     split(.$endpoint)
#   outd <- run_curvep_batch(x[[1]],
#                          directionality = 1,
#                          n_sample = 2,
#                          threshold = seq(5, 10, by = 5),
#                          other_paras = list(CARR = 20, TrustHi = TRUE))
#
#   outd <- outd %>% dplyr::select(threshold, direction, repeat_id) %>% dplyr::distinct()
#   expect_true(nrow(outd) == 4, info = "a failed run" )})
#
#
# test_that("resp boot error dataset job run", {
#
#   error <- rnorm(100, 0, 10)
#   data("zfishbeh")
#   x <- zfishbeh %>%
#     split(.$endpoint)
#   outd <- run_curvep_batch(x[[1]],
#                          directionality = 1,
#                          n_sample = 2,
#                          threshold = seq(5, 10, by = 5),
#                          other_paras = list(CARR = 20, TrustHi = TRUE), vehicle_data = error )
#
#   outd <- outd %>% dplyr::select(threshold, direction, repeat_id) %>% dplyr::distinct()
#   expect_true(nrow(outd) == 4, info = "a failed run" )
# })
