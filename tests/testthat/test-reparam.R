context("reparameters")

test_that("change CARR", {
  data("zfishbeh")
  x <- zfishbeh %>%
    split(.$endpoint)
  outd <- run_curvep_job(x[[1]],
                         directionality = 1,
                         n_sample = NULL,
                         threshold = 15,
                         other_paras = list(CARR = 20, TrustHi = TRUE))
  other_paras <- list(CARR = 50, TrustHi = TRUE)
  outd2 <- reparam_curvep_job(outd, other_paras = other_paras)
  expect_true(outd2[1,]$input[[1]]$paras[['CARR']] == 50)
})

test_that("change directionality", {
  data("zfishbeh")
  x <- zfishbeh %>%
    split(.$endpoint)
  outd <- run_curvep_job(x[[1]],
                         directionality = 1,
                         n_sample = NULL,
                         threshold = 15,
                         other_paras = list(CARR = 20, TrustHi = TRUE))
  outd2 <- reparam_curvep_job(outd, directionality = -1)
  expect_true(outd2[1,]$input[[1]]$paras[['RNGE']] == -1e+06)
})


test_that("wrong directionality", {
  data("zfishbeh")
  x <- zfishbeh %>%
    split(.$endpoint)
  outd <- run_curvep_job(x[[1]],
                         directionality = 1,
                         n_sample = NULL,
                         threshold = 15,
                         other_paras = list(CARR = 20, TrustHi = TRUE))
  expect_error(
    outd2 <- reparam_curvep_job(outd, directionality = c(-1,1))
  )

})

