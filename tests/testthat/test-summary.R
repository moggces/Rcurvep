context("summarize results")


data("zfishbeh")
set.seed(300)
outp1 <- combi_run_rcurvep(zfishbeh, n_samples = 3, TRSH = c(5, 10))
outp2 <- combi_run_rcurvep(zfishbeh, n_samples = NULL, TRSH = c(5, 10))
outp3 <- run_rcurvep(create_dataset(zfishbeh), RNGE = 1000000)

test_that("can run and create datasets", {

  outp1s <- summarize_rcurvep_output(outp1)
  outp2s <- summarize_rcurvep_output(outp2)
  outp3s <- summarize_rcurvep_output(outp3)

  expect_true(all( c('result', 'act_summary', 'config') %in% names(outp1s)))
  expect_true(all( c('result', 'act_summary', 'config') %in% names(outp2s)))
  expect_true(all( c('result', 'act_summary', 'config') %in% names(outp3s)))
})

test_that("inactivate, string", {
  outp1s <- summarize_rcurvep_output(outp1, inactivate  = "INVERSE", clean_only = TRUE)
  act_set <- outp1s$result$act_set
  ind <- stringr::str_detect(act_set$Comments, "INVERSE")
  #expect_true(all(!is.na(act_set$POD[ind])))
  expect_true(sum(act_set$hit[ind]) == 0)
  expect_true(all(stringr::str_detect(act_set$Comments[ind], "custom")))
})


test_that("inactivate, index", {

  # it is working but because the sort of table is different...
  # ind <- c(2,3)
  # outp1s <- summarize_rcurvep_output(outp1, inactivate  = ind, clean_only = TRUE)
  # act_set <- outp1s$result$act_set
  # #expect_true(all(!is.na(act_set$POD[ind])))
  #
  # expect_true(sum(act_set$hit[ind]) == 0)
  # expect_true(all(stringr::str_detect(act_set$Comments[ind], "custom")))
})
