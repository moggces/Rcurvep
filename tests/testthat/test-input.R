context("input parameters")


test_that("threshold is not numeric", {
  data("zfishbeh")
  x <- zfishbeh %>%
    dplyr::group_by(endpoint, chemical, concs) %>%
    dplyr::slice(1) %>%
    split(.$endpoint)
  expect_error(
    run_curvep_job(x[[1]],
                         directionality = 0,
                         n_sample = NULL,
                         threshold = "15",
                         other_paras = list())
  )
})


test_that("threshold is not numeric in a list", {
  data("zfishbeh")
  x <- zfishbeh %>%
    dplyr::group_by(endpoint, chemical, concs) %>%
    dplyr::slice(1) %>%
    split(.$endpoint)
  expect_error(
    run_curvep_job(x[[1]],
                   directionality = 0,
                   n_sample = NULL,
                   threshold = list("-1" = c(5), "1" = c("a", "b")),
                   other_paras = list())
  )
})


test_that("threshold and direction do not match", {
  data("zfishbeh")
  x <- zfishbeh %>%
    dplyr::group_by(endpoint, chemical, concs) %>%
    dplyr::slice(1) %>%
    split(.$endpoint)
  expect_error(
    run_curvep_job(x[[1]],
                   directionality = 0,
                   n_sample = NULL,
                   threshold = list("-1" = c(5), "a" = c(5, 10)),
                   other_paras = list())
  )
})

