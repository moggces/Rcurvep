context("BMR")

data(zfishdev_act)
sumd <- summarize_rcurvep_output(zfishdev_act, clean_only = TRUE)

test_that("can run", {

  bmr_out <- estimate_dataset_bmr(sumd, plot = FALSE)
  expect_true(all(names(bmr_out) %in% c('stats', 'outcome')))
})

test_that("set p1 and p2", {
  bmr_out <- estimate_dataset_bmr(sumd, p1 = 2, p2 = 18, plot = FALSE)
  statsd <- bmr_out$stats
  expect_true(all(unique(statsd$p1_ori) >= 2, unique(statsd$p1_exp) == 2,
                  unique(statsd$p2_ori) == 18, unique(statsd$p2_exp) == 18))
})

