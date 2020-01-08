context("input parameters")

test_that("check parameter names", {

  defaults <- curvep_defaults()

  # default
  config <- .check_config_name(config = defaults, RNGE = 1000000)
  expect_equal(defaults$MXDV, config$MXDV)

  # add through command
  config <- .check_config_name(config = defaults, MXDV = 10)
  expect_equal(config$MXDV, 10)

  # add through config
  custom <- defaults
  custom[['MXDV']] <- 10
  config <- .check_config_name(config = custom)
  expect_equal(config$MXDV, 10)

  # add non existent parameter
  custom[['NON']] <- 10
  expect_error(.check_config_name(config = custom))

  # try to overwrite RNGE
  config <- .check_config_name(config = defaults, RNGE = 100, RNGE = 100000)
  expect_equal(config$RNGE, 100)

})

dats <- zfishbeh %>%
  dplyr::group_by(endpoint, chemical, conc) %>%
  dplyr::slice(1) %>%
  dplyr::ungroup() %>%
  split(.$endpoint)

test_that("check mask value out of range", {

  expect_error(.check_mask_input(c(1, 6), dats[[1]]))
})

test_that("mask column available", {
  dats1 <- dats[[1]]
  dats1$mask <- 0
  expect_warning(.check_mask_input(c(1, 6), dats1))
})


test_that("check parameter values", {

  defaults <- curvep_defaults()

  config <- .check_config_name(config = defaults, RNGE = "1000000")
  expect_error(.check_config_value(config))

})

test_that("check input basic dataset", {

  data("zfishbeh")

  # too many resps per conc
  expect_error(.check_dat_base(zfishbeh))

  # one resp per conc
  dats <- zfishbeh %>%
    dplyr::group_by(endpoint, chemical, conc) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    split(.$endpoint)
  expect_equal(.check_dat_base(dats[[1]]), dats[[1]])

  # mask
  d <- dats[[1]] %>% dplyr::mutate(mask = 0)
  expect_equal(.check_dat_base(d), d)

  # extra unneeded columns
  d <- dats[[1]] %>% dplyr::mutate(sample_id = 0)
  expect_equal(.check_dat_base(d), dats[[1]])


})

test_that("input two types of dataset", {
  expect_type(.check_dat(tibble::tibble(endpoint = 1, chemical = 1, resp = 1, conc = 1, mask = 1)), "list")
  expect_error(.check_dat(tibble::tibble(endpoint = 1, chemical = 1, resp = 1, mask = 1)))
  expect_type(.check_dat(tibble::tibble(endpoint = 1, chemical = 1, n_in = 1, N = 1, conc = 1, mask = 1)), "list")
})

test_that("vehicle control data",  {
  expect_type(.check_vdata(rnorm(10), "continuous"), "double")
  expect_error(.check_vdata(c(1, 10, 100, 50, NA), "continuous"))
  expect_error(.check_vdata(rnorm(10), "dichotomous"))
})

test_that("keep_sets argument", {
  expect_error(.check_keep_sets(c("act_set", "xx"), c("act_set", "resp_set", "fp_set"), "act_set"))
  expect_length(.check_keep_sets(c("act_set", "act_set", "fp_set"),  c("act_set", "resp_set", "fp_set"), "act_set"), 2)
})

test_that("check_result_sets", {
  expect_error(.check_result_sets(c("act_set")))
  outp <- run_rcurvep(create_dataset(zfishbeh))
  expect_length(.check_result_sets(outp$result), 3)
  expect_error(.check_result_sets(outp$result['fp_set']))
})

