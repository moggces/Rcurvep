context("BMR")

test_that("BMR selection", {
  data("zfishdev_act")
  zfish2 <-
    tidyr::separate(zfishdev_act, .data$dduid, c("endpoint", "chemical", "directionality"), sep = "#")
  zfish2 <- dplyr::mutate(zfish2, POD = ifelse(is.na(POD), conc_highest, POD))
  result <- identify_basenoise_threshold(zfish2, id = c("endpoint", "directionality"),
       chemical = "chemical", threshold = "thres", potency = "POD" )
  expect_true(tibble::is_tibble(result))

})
