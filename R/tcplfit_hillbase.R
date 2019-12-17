
fit_modls <- function(Conc, Resp, Mask = NULL, modls = c("hill", "cnst"), ...) {
  args <- list(...)
  result <- purrr::map(
    modls, ~ do.call(
      fit_modl_in, c(list(Conc = Conc, Resp = Resp, Mask = Mask, modl = .x), args)
    )
  ) %>% rlang::set_names(modls)
  result <- sort_modl_aic(result)
  return(result)
}

fit_modl_in <- function(Conc, Resp, Mask, modl, ...) {
  l <- list(...)
  func_name <- stringr::str_c("fit_", modl, "_modl")
  pref <- stringr::str_glue("^{modl}_")
  l_p <- l[stringr::str_detect(names(l), pref)]
  l_p <- l_p %>% rlang::set_names(stringr::str_remove(names(l_p), pref))

  result <- NULL
  if (modl == "cnst") {
    result <- do.call(
      get(eval(func_name)), c(list(Resp = Resp), l_p)
    )
  } else if (modl == "hill") {
    result <- do.call(
      get(eval(func_name)), c(list(Conc = Conc, Resp = Resp, Mask = Mask), l_p)
    )
  }

  return(result)
}


#' Fit Hill model
#'
#' @param Conc
#' @param Resp
#' @param Mask
#' @param pdir preferred direction, 1 (increasing) or -1 (decreasing), or both then AIC is used to select one
#' @param ... can be used to replace the initial config of Hill fit (e.g., f )
#'
#' @return a (selected) Hill model fit, a list including modl, fit, ga, gw, tp
#' tp < 0 means that decreasing
#' @export

fit_hill_modl <- function(Conc, Resp, Mask = NULL, pdir = c(1, -1), ...) {

  # need to test the input ....
  if (length(pdir) > 1) {
    result <- purrr::map(
      pdir, ~ fit_hill_modl_in(Conc = Conc, Resp = Resp, Mask = Mask, pdir = .x, ...))
  } else {
    result <- list(
      fit_hill_modl_in(Conc = Conc, Resp = Resp, Mask = Mask, pdir = pdir, ...))
  }
  result <- select_modl_aic(result)
  return(result)
}

select_modl_aic <- function(lmodls) {

  if(length(lmodls) > 1) {
    aics <- purrr::map_dbl(lmodls, ~ .x[['aic']])
    result <- lmodls[[which.min(aics)]]
  } else {
    result <- lmodls[[1]]
  }
  return(result)
}

sort_modl_aic <- function(lmodls) {

    aics <- purrr::map_dbl(lmodls, ~ .x[['aic']])
    result <- lmodls[order(aics)]

  return(result)
}

#' Fit Hill
#'
#' adjust the direction
#'
#' @param Conc
#' @param Resp
#' @param Mask
#' @param pdir allow only one
#' @param ...
#'
#' @return
#' @keywords internal
#'
#'
fit_hill_modl_in <- function(Conc, Resp, Mask = NULL, pdir = c(1, -1), ...) {

  args <- list(...)

  # default values
  hill <- 0L
  haic <- NA_real_
  hill_tp <- NA_real_
  hill_ga <- NA_real_
  hill_gw <- NA_real_
  hill_er <- NA_real_


  # need to check Resp and Mask has the same length
  cleand <- mask_resp_conc(Conc = Conc, Resp = Resp, Mask = Mask)
  conc <- cleand$Conc
  resp <- cleand$Resp

  # adjust the direction
  if (pdir == -1) resp <- resp*-1

  # replace the default config if needed
  config <- get_hill_fit_config(conc, resp)
  config <- modifyList(config, args)


  # Input for optimization
  h <- config$theta
  hUi <- config$ui
  hCi <- config$ci
  f <- config$f

  # optimize the hill model
  hfit <- try(constrOptim(theta = h,
                          f = get(eval(f)),
                          ui = hUi,
                          ci = hCi,
                          mu = 1e-6,
                          method = "Nelder-Mead",
                          control = list(fnscale = -1,
                                         reltol = 1e-10,
                                         maxit = 6000),
                          lconc = conc,
                          resp = resp),
              silent = TRUE)

  # make sure there is a fit
  if (!is(hfit, "try-error")) { # Hill model fit the data
    hill <- 1L
    haic <- 8 - 2*hfit$value # 2*length(hfit$par) - 2*hfit$value

    hill_tp <- hfit$par[1]
    hill_ga <- hfit$par[2]
    hill_gw <- hfit$par[3]
    hill_er <- hfit$par[4]

    # adjust the direction
    if (pdir == -1) hill_tp <- hill_tp*-1

  }

  # output
  result <- list(
    modl = "hill",
    fit = hill,
    aic = haic,
    tp = hill_tp,
    ga = hill_ga,
    gw = hill_gw,
    er = hill_er
  )
  return(result)
}



#' Get the default config for the hill fit
#'
#' @param Conc a vector of log10 conc
#' @param Resp a vector of responses
#' @param optimf default is tcplObjHill but can change to ObjHillnorm
#'
#' @return a list of theta (inital values), f (optimf), ui (bound matrix), ci (bounds)
#'
#' @export
#'
get_hill_fit_config <- function(Conc, Resp, optimf = "tcplObjHill" ) {

  # used values
  rmds <- tapply(Resp, Conc, median)
  mmed <- max(rmds)
  mmed_conc <- as.numeric(names(which.max(rmds)))
  resp_max <- max(Resp)
  logc_min <- min(Conc)
  logc_max <- max(Conc)
  LBratio <- 1/100 # will make the lower bound of ga even lower (log10(1/100))

  # scale term estimate
  er_est <- get_est_error(Resp)

  # contraint the estimate (need to compare with constraint matrix)
  hbnds <- c(0, -1.2*resp_max, # tp bounds
             logc_min +log10(LBratio), -(logc_max + 0.5), # ga bounds
             0.3, -8) # gw bounds

  # starting parameters
  h <- c(mmed, # top (tp)
         mmed_conc - 1, # logAC50 (ga)
         1.2, # hill coefficient (gw)
         er_est) # logSigma (er)

  # adjust some starting parameters
  if (h[1] == 0) h[1] <- 0.1
  # adjust guess for k to be within bounds
  if (h[2] <= hbnds[3]) {
    h[2] <- hbnds[3] + max(0.0001*abs(hbnds[3]),0.000001) #2nd part is to handle bounds of log10(1)=0
  }

  # constraint matrix (k x p)
  #                tp   ga   gw   er
  hUi <- matrix(c( 1,   0,   0,   0,
                   -1,   0,   0,   0,
                   0,   1,   0,   0,
                   0,  -1,   0,   0,
                   0,   0,   1,   0,
                   0,   0,  -1,   0),
                byrow = TRUE, nrow = 6, ncol = 4)

  # constraint vector of length k
  hCi <- matrix(hbnds, nrow = 6, ncol = 1)


  # output
  result <- list(
    theta = h,
    f = optimf,
    ui = hUi,
    ci = hCi
  )
  return(result)
}



#' Guess the initial scaling factor
#'
#' @param Resp a vector of resps
#'
#' @return a initial scaling factor
#' @keywords internal

get_est_error <- function(Resp) {
  rmad <- mad(Resp)
  er_est <- if (rmad > 0) log(rmad) else log(1e-32)
  return(er_est)
}



#' Adjust the conc and resp using the mask
#'
#' @param Conc a vector of log conc
#' @param Resp a vector of resp
#' @param Mask a vector of 1/0 for masking
#'
#' @return  a list of updated Conc and Resp
#' @keywords internal
#'
mask_resp_conc <- function(Conc, Resp, Mask) {

  if (!is.null(Mask)) {
    Resp <- replace(Resp, as.logical(Mask), NA)
    ind <- which(!is.na(Resp))
    Resp <- Resp[ind]
    Conc <- Conc[ind]
  }

  return(list(Conc = Conc, Resp = Resp))
}




#' Fit a constant model using object function tcplObjCnst
#'
#' @param Resp a vector of responses
#' @param optimf object function, tcplObjCnst
#'
#' @return a list of parameters, modl, fit (1,0), aic, er
#' @keywords internal
#'
#'
fit_cnst_modl <- function(Resp, optimf = "tcplObjCnst") {

  # default value
  cnst <- 0L
  cnst_er <- NA_real_
  caic <- NA_real_


  # scale term estimate
  er_est <- get_est_error(Resp)

  #optimization
  cfit <- optim(er_est,
                fn = get(eval(optimf)),
                method = "Brent",
                lower = er_est - 2,
                upper = er_est + 2,
                control = list(fnscale = -1,
                               reltol = 1e-4,
                               maxit = 500),
                resp = Resp)

  # check if there is a fit
  if (!is(cfit, "try-error")) {
    cnst <- 1L
    cnst_er <- cfit$par
    caic <- 2 - 2*cfit$value # 2*length(cfit$par) - 2*cfit$value
  }

  # output
  result <- list(
    modl = "cnst",
    fit = cnst,
    er = cnst_er,
    aic = caic
  )

  return(result)

}
