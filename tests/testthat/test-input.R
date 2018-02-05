context("input parameters")


test_that("no para input", {
  data("zfishbeh")
  x <- zfishbeh %>%
    group_by(endpoint, chemical, concs) %>%
    slice(1) %>%
    split(.$endpoint)
  outd <- run_curvep_job(x[[1]],
                         directionality = 0,
                         n_sample = NULL,
                         threshold = 15,
                         other_paras = list())
  expect_true(sum(colnames(outd) %in% c('input', 'output', 'activity')) == 3, info = "a failed run" )
})
