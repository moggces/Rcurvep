context("BMR")

test_that("NA in potency column", {
  data("zfishdev_act")

  expect_error(
    identify_basenoise_threshold(zfishdev_act), "NA is not allowed in the potency column",
    ignore.case = TRUE
  )

})

test_that("BMR selection", {
  data("zfishdev_act")
  zfish2 <- dplyr::mutate(zfishdev_act, POD = ifelse(is.na(POD), conc_highest, POD))
  result <- identify_basenoise_threshold(zfish2, endpoint = "endpoint", direction = "direction",
       chemical = "chemical", threshold = "threshold", potency = "POD" )
  expect_true(tibble::is_tibble(result))

})

test_that("user selected p2", {
  data("zfishdev_act")
  zfish2 <- dplyr::mutate(zfishdev_act, POD = ifelse(is.na(POD), conc_highest, POD))
  result <- identify_basenoise_threshold(
    zfish2, endpoint = "endpoint", direction = "direction",
            chemical = "chemical", threshold = "threshold", potency = "POD", p2 = 7)
  expect_true(unique(result$p2) == 7)

})


test_that("thresDist flag", {
  data("zfishdev_act")
  zfish2 <- dplyr::mutate(zfishdev_act, POD = ifelse(is.na(POD), conc_highest, POD))
  result <- identify_basenoise_threshold(
    zfish2, endpoint = "endpoint", direction = "direction",
    chemical = "chemical", threshold = "threshold", potency = "POD")
  comment <- dplyr::filter(result, thresDist == 1) %>% dplyr::pull(thresDistComment)

  expect_true(sum(comment %in% c("OK", "inadvisable")) == 2)

})
