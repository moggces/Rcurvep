context("output")

test_that("output activities", {
  data("zfishdev")
  x <- zfishdev %>%
    split(.$endpoint)
  outd <- run_curvep_job(x[[1]],
                         directionality = 1,
                         n_sample = 1,
                         threshold = 15,
                         other_paras = list(CARR = 20, TrustHi = TRUE))
  acts <- extract_curvep_data(outd, "act")
  expect_true(sum(colnames(acts) %in% c('POD', 'wAUC', 'Emax')) == 3, info = "a failed run" )

})

test_that("id columns are in the output", {
  data("zfishbeh")
  x <- zfishbeh %>%
    split(.$endpoint)
  outd <- run_curvep_job(x[[1]],
                         directionality = 0,
                         n_sample = NULL,
                         threshold = list("1" = c(15,20), "-1" = 30),
                         other_paras = list(CARR = 20, TrustHi = TRUE))
  expect_true(sum(colnames(outd) %in% c('endpoint', 'chemical', 'direction', 'threshold')) == 4, info = "a failed run" )
})
