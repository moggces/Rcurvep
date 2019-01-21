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


test_that("output activities summary", {
  data("zfishdev")
  x <- zfishdev %>%
    split(.$endpoint)
  outd <- run_curvep_job(x[[1]],
                         directionality = 1,
                         n_sample = 10,
                         threshold = c(15, 20),
                         other_paras = list(CARR = 20, TrustHi = TRUE))
  summ <- extract_curvep_data(outd, "summary")
  expect_true(sum(colnames(summ) %in% c('POD_med', 'POD_ciu', 'POD_cil', 'hit_confidence', "concs", "resps", "resps_in")) == 7, info = "a failed run" )

})

test_that("output activities summary (no simulation)", {
  data("zfishdev")
  x <- zfishdev %>%
    split(.$endpoint)
  outd <- run_curvep_job(x[[1]],
                         directionality = 1,
                         n_sample = NULL,
                         threshold = c(15, 20),
                         other_paras = list(CARR = 20, TrustHi = TRUE))
  summ <- extract_curvep_data(outd, "summary")
  expect_true(sum(colnames(summ) %in% c('POD_med', 'POD_ciu', 'POD_cil', 'hit_confidence', "concs", "resps", "resps_in")) == 7, info = "a failed run" )

})


test_that("remove comment activities", {
  data("zfishdev")
  x <- zfishdev %>%
    split(.$endpoint)
  outd <- run_curvep_job(x[[1]],
                         directionality = 1,
                         n_sample = 10,
                         threshold = 5,
                         other_paras = list(CARR = 20, TrustHi = TRUE))
  acts <- extract_curvep_data(outd, "act", modifier = "CHECK")
  expect_true(sum(acts$hit[grepl("CHECK",acts$Comments)]) == 0)

})

test_that("summarize curvep act data", {
  data("zfishdev")
  x <- zfishdev %>%
    split(.$endpoint)
  outd <- run_curvep_job(x[[1]],
                         directionality = 1,
                         n_sample = 10,
                         threshold = 5,
                         other_paras = list(CARR = 20, TrustHi = TRUE), simplify_output = TRUE)
  #outd <- outd %>% select(-wAUC_prev) #it will not work
  acts <- summarize_curvep_output(outd, col_names = c("POD", "EC50"), modifier = "INVERSE", conf_level = 0.9)
  expect_true(ncol(acts) == 11)
})
