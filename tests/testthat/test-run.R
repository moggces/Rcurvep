context("running job")


test_that("percent data with n_sample = NULL", {
  data("zfishdev")
  x <- zfishdev %>%
    split(.$endpoint)
  outd <- run_curvep_job(x[[1]],
                         directionality = 1,
                         n_sample = NULL,
                         threshold = 15,
                         other_paras = list(CARR = 20, TrustHi = TRUE))
  expect_true(sum(colnames(outd) %in% c('input', 'output', 'activity')) == 3, info = "a failed run" )
})


test_that("resp data with n_sample = NULL", {
  data("zfishbeh")
  x <- zfishbeh %>%
    split(.$endpoint)
  outd <- run_curvep_job(x[[1]],
                         directionality = 1,
                         n_sample = NULL,
                         threshold = 15,
                         other_paras = list(CARR = 20, TrustHi = TRUE))
  expect_true(sum(colnames(outd) %in% c('input', 'output', 'activity')) == 3, info = "a failed run" )
})


test_that("resp data original with multiple thresholds", {
  data("zfishbeh")
  x <- zfishbeh %>%
    split(.$endpoint)
  outd <- run_curvep_job(x[[1]],
                         directionality = 0,
                         n_sample = NULL,
                         threshold = list("1" = c(15,20), "-1" = c(30)),
                         other_paras = list(CARR = 20, TrustHi = TRUE))

  outd <- outd %>% dplyr::select(threshold, direction) %>% dplyr::distinct() %>%
    dplyr::arrange(direction, threshold) %>% dplyr::pull(threshold)
  expect_equal(outd, c(30, 15, 20))
})


test_that("percent boot dataset job run", {
  data("zfishdev")
  x <- zfishdev %>%
    split(.$endpoint)
  outd <- run_curvep_job(x[[1]],
                         directionality = 1,
                         n_sample = 2,
                         threshold = seq(5, 10, by = 5),
                         other_paras = list(CARR = 20, TrustHi = TRUE))

  outd <- outd %>% dplyr::select(threshold, direction, repeat_id) %>% dplyr::distinct()

  expect_true(nrow(outd) == 4, info = "a failed run" )
})


test_that("resp boot dataset job run", {
  data("zfishbeh")
  x <- zfishbeh %>%
    split(.$endpoint)
  outd <- run_curvep_job(x[[1]],
                         directionality = 1,
                         n_sample = 2,
                         threshold = seq(5, 10, by = 5),
                         other_paras = list(CARR = 20, TrustHi = TRUE))

  outd <- outd %>% dplyr::select(threshold, direction, repeat_id) %>% dplyr::distinct()
  expect_true(nrow(outd) == 4, info = "a failed run" )})


test_that("resp boot error dataset job run", {

  error <- rnorm(100, 0, 10)
  data("zfishbeh")
  x <- zfishbeh %>%
    split(.$endpoint)
  outd <- run_curvep_job(x[[1]],
                         directionality = 1,
                         n_sample = 2,
                         threshold = seq(5, 10, by = 5),
                         other_paras = list(CARR = 20, TrustHi = TRUE), vehicle_data = error )

  outd <- outd %>% dplyr::select(threshold, direction, repeat_id) %>% dplyr::distinct()
  expect_true(nrow(outd) == 4, info = "a failed run" )
})
