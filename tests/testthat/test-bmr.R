context("BMR")

test_that("NA in potency column", {
  data("zfishdev_act")
  zfish2 <- zfishdev_act %>% dplyr::select(-conc_highest)
  expect_error(
    get_baseline_threshold(zfish2), "NA is not allowed in the potency column",
    ignore.case = TRUE
  )

})

test_that("run", {
  data("zfishdev_act")
  result <- get_baseline_threshold(zfishdev_act)
  expect_true(rlang::is_list(result))

})

test_that("run with non default column names", {
  data("zfishdev_act")
  zfish2 <- zfishdev_act %>%
    dplyr::rename(BMC = POD)
  result <- get_baseline_threshold(zfish2, potency = "BMC")
  expect_true(rlang::is_list(result))

})

test_that("user selected p2", {
  data("zfishdev_act")
  zfish2 <- zfishdev_act %>% dplyr::filter(endpoint == "percent_affected_96")
  result <- get_baseline_threshold(zfish2)
  result <- cal_knee_point(result[[1]], "threshold", "pooled_variance", p2_raw = 7)
  expect_true(unique(result[[2]]$p2_raw) == 7)

})


test_that("thresDist flag", {
  data("zfishdev_act")
  result <- get_baseline_threshold(zfishdev_act)
  comment <- result[[2]] %>% dplyr::pull(thresComment)

  expect_true(sum(comment %in% c("OK", "check")) == 2)

})
