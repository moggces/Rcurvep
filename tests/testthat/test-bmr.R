context("BMR")

test_that("NA in potency column", {
  data("zfishdev_act")

  expect_error(
    identify_basenoise_threshold(zfishdev_act), "NA is not allowed in the potency column",
    ignore.case = TRUE
  )

})

test_that("run", {
  data("zfishdev_act")
  zfish2 <- dplyr::mutate(zfishdev_act, POD = ifelse(is.na(POD), conc_highest, POD))
  result <- identify_basenoise_threshold(zfish2 )
  expect_true(rlang::is_list(result))

})

test_that("run with non default column names", {
  data("zfishdev_act")
  zfish2 <- dplyr::mutate(zfishdev_act, POD = ifelse(is.na(POD), conc_highest, POD)) %>%
    dplyr::rename(BMC = POD)
  result <- identify_basenoise_threshold(zfish2, potency = "BMC" )
  expect_true(rlang::is_list(result))

})

test_that("user selected p2", {
  data("zfishdev_act")
  zfish2 <- dplyr::mutate(zfishdev_act, POD = ifelse(is.na(POD), conc_highest, POD))
  zfish2 <- zfish2 %>% dplyr::filter(endpoint == "percent_affected_96")
  result <- identify_basenoise_threshold(zfish2)
  result <- cal_knee_point(result[[1]], "threshold", "pooled_variance", p2_raw = 7)
  expect_true(unique(result[[2]]$p2_raw) == 7)

})


test_that("thresDist flag", {
  data("zfishdev_act")
  zfish2 <- dplyr::mutate(zfishdev_act, POD = ifelse(is.na(POD), conc_highest, POD))
  result <- identify_basenoise_threshold(zfish2)
  comment <- result[[2]] %>% dplyr::pull(thresComment)

  expect_true(sum(comment %in% c("OK", "check")) == 2)

})
