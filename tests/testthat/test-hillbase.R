context("model fit using one chemical")

concd <- c(-9, -8, -7, -6, -5, -4)
respd_up <- c(0, 2, 30, 40, 50, 60)
respd_down <- c(0, 2, 30, 40, 50, 60)*-1


test_that("use one model", {

  outp <- fit_modls(concd, respd_up, modls = c("cnst"))
  expect_true(length(outp) == 1 && outp[[1]]$modl == "cnst")

})

test_that("use hill + one direction", {
  outp <- fit_modls(concd, respd_up, hill_pdir = -1, modls = "hill")
  expect_true(is.na(outp$hill$tp))
})

test_that("use hill + a different object function", {
  outp1 <- fit_modls(concd, respd_up, modls = "hill", hill_f = "ObjHillnorm")
  outp2 <- fit_modls(concd, respd_up, modls = "hill")

  expect_true(outp1$hill$ga != outp2$hill$ga)
})

test_that("use cc2", {
  outp <- fit_modls(concd, respd_up, modls = "cc2")
  expect_true(outp$cc2$tp > 0)
})

test_that("use cc2 + change of classSD", {
  outp <- fit_modls(concd, respd_up, modls = "cc2", cc2_classSD = 30)
  expect_true(outp$cc2$cc2 == 4)
})


test_that("warnings", {
  expect_warning(fit_modls(concd, respd_up, modls = "cnst", hill_f = "ObjHillnorm"))
  expect_warning(fit_modls(concd, respd_up, modls = "hill", hill_fx = "ObjHillnorm"))
  expect_warning(fit_modls(concd, respd_up, modls = "hill", xx = "ss"))

})


test_that("errors", {
  expect_error(fit_modls(concd, respd_up[-1]), modls = "hill")
  expect_error(fit_modls(concd[1:3], respd_up[1:3]), modls = "cc2")
  expect_error(fit_modls(concd, respd_up, Mask = rep(1, length(concd)),  modls = "hill"))
  expect_error(fit_modls(c(concd, -3), c(respd_up, NA), modls = "hill"))

})



