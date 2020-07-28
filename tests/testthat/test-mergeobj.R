context("merge rcurvep objects")

data("zfishbeh")
dat <- zfishbeh

# decreasing/no sampling
out1 <- combi_run_rcurvep(
  dat, n_samples = NULL, TRSH = 10,  keep_sets = c("act_set", "resp_set"))

# decreasing/no sampling + mask
out2 <- combi_run_rcurvep(
  dat, n_samples = NULL, TRSH = 10, mask = 1,  keep_sets =  c("act_set", "resp_set"))

# increasing/no sampling
out3 <- combi_run_rcurvep(
  dat, n_samples = NULL, TRSH = 10, mask = 1, RNGE = 1000000, keep_sets =  c("act_set", "resp_set"))

# decreasing/sampling
out4 <- combi_run_rcurvep(
  dat, n_samples = 2, TRSH = 10, seed = 300)

# increasing/sampling
out5 <- combi_run_rcurvep(
  dat, n_samples = 2, TRSH = 10, RNGE = 100000, seed = 300)

# hill fit
out6 <- run_fit(zfishbeh, keep_sets = "fit_set")

# decreasing/sampling but only act_set
set.seed(300)
out7 <- combi_run_rcurvep(
  dat, n_samples = 2, TRSH = 10, keep_sets = c("act_set", "resp_set"))

test_that("cannot merge + weird merge", {

  # curvpe != hill
  expect_error(merge_rcurvep_obj(out1, out6))
  # size of lsets
  expect_error(merge_rcurvep_obj(out1, out4))
  # can merge but may not meaningful
  expect_length(merge_rcurvep_obj(out1, out7), 2)
})


test_that("merge objects", {

  m1 <- merge_rcurvep_obj(out1, out2, out3)
  expect_true(nrow(m1$result$act_set) == 22)


  m2 <- merge_rcurvep_obj(out4, out5)
  expect_true(nrow(m2$result$act_set) == 44)

})

test_that("summarize the merged object", {


  # can merge
  m1 <- merge_rcurvep_obj(out1, out2, out3)
  out <- summarize_rcurvep_output(m1, inactivate = c(2, 3), clean_only = TRUE)
  expect_length(out, 2)

  # can merge
  m2 <- merge_rcurvep_obj(out4,out5)
  out <- summarize_rcurvep_output(m2, inactivate = c(2, 3), clean_only = TRUE)
  expect_length(out, 2)

  ## note the results is not 22 anymore because the base cols are different
  out <- summarize_rcurvep_output(m2, inactivate = c(2, 3), clean_only = TRUE)
  expect_true(nrow(m2$result$act_set) != 22)
})

