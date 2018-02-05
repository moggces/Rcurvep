context("running job")

test_that("resp original dataset job run", {
  data("zfishbeh")
  x <- zfishbeh %>%
    group_by(endpoint, chemical, concs) %>%
    slice(1) %>%
    split(.$endpoint)
  outd <- run_curvep_job(x[[1]],
                         directionality = 0,
                         n_sample = NULL,
                         threshold = 15,
                         other_paras = list(CARR = 20, TrustHi = TRUE))
  expect_true(sum(colnames(outd) %in% c('input', 'output', 'activity')) == 3, info = "a failed run" )
})


test_that("resp original dataset median job run", {
  data("zfishbeh")
  x <- zfishbeh %>%
    split(.$endpoint)
  outd <- run_curvep_job(x[[1]],
                         directionality = 0,
                         n_sample = NULL,
                         threshold = 15,
                         other_paras = list(CARR = 20, TrustHi = TRUE))
  expect_true(sum(colnames(outd) %in% c('input', 'output', 'activity')) == 3, info = "a failed run" )
})

test_that("percent boot dataset job run", {
  data("zfishdev")
  x <- zfishdev %>%
    split(.$endpoint)
  outd <- run_curvep_job(x[[1]],
                         directionality = 1,
                         n_sample = 1,
                         threshold = seq(5, 10, by = 5),
                         other_paras = list(CARR = 20, TrustHi = TRUE))
  expect_true(sum(colnames(outd) %in% c('input', 'output', 'activity')) == 3, info = "a failed run" )
})


test_that("resp boot dataset job run", {
  data("zfishbeh")
  x <- zfishbeh %>%
    split(.$endpoint)
  outd <- run_curvep_job(x[[1]],
                         directionality = 0,
                         n_sample = 2,
                         threshold = seq(5, 10, by = 5),
                         other_paras = list(CARR = 20, TrustHi = TRUE))
  expect_true(sum(colnames(outd) %in% c('input', 'output', 'activity')) == 3, info = "a failed run" )
})


test_that("resp boot error dataset job run", {

  error <- rnorm(100, 0, 10)
  data("zfishbeh")
  x <- zfishbeh %>%
    split(.$endpoint)
  outd <- run_curvep_job(x[[1]],
                         directionality = 0,
                         n_sample = 2,
                         threshold = seq(5, 10, by = 5),
                         other_paras = list(CARR = 20, TrustHi = TRUE), vehicle_data = error )
  expect_true(sum(colnames(outd) %in% c('input', 'output', 'activity')) == 3, info = "a failed run" )
})
