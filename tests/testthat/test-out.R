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
