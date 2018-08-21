context("BMR")

test_that("BMR selection", {
  data("zfishdev_act")
  zfish2 <- dplyr::mutate(zfishdev_act, POD = ifelse(is.na(POD), conc_highest, POD))
  result <- identify_basenoise_threshold(zfish2, endpoint = "endpoint", direction = "direction",
       chemical = "chemical", threshold = "threshold", potency = "POD" )
  expect_true(tibble::is_tibble(result))

})

test_that("BMR selection", {
  data("zfishdev_act")
  zfish2 <- dplyr::mutate(zfishdev_act, POD = ifelse(is.na(POD), conc_highest, POD))
  result <- identify_basenoise_threshold(
    zfish2, endpoint = "endpoint", direction = "direction",
            chemical = "chemical", threshold = "threshold", potency = "POD", p2 = 7)
  expect_true(unique(result$p2) == 7)

})
