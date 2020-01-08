
#' Fit one concentration-response data using types of models
#'
#' A convenient function to fit data using available models
#' and to sort the outcomes by AIC values.
#'
#'
#'
#' @details
#' The backbone of fit using hill and cnst is based on the implementation from tcpl package.
#' But the lower bound of ga is lower by log10(1/100).
#'
#'

#' @param Conc A vector of log10 concentrations.
#' @param Resp A vector of numeric responses.
#' @param Mask Default = NULL or a vector of 1 or 0.
#'   1 is for masking the respective response.
#' @param modls The model types for the fitting. Multiple values are allowed.
#'   Currently Hill model (hill) and constant model (cnst) are implemented.
#'   Default = c("hill", "cnst").
#' @param ... The named input configurations for replacing the default configurations.
#'   The input configuration needs to add model type as the prefix.
#'   For example, hill_pdir = -1 will set the Hill fit only to the decreasing direction.
#'
#' @return A list of components named by the models.
#'   The models are sorted by their AIC values. Thus, the first component has the best fit.
#'
#'  ## hill
#'   Fit output from Hill equation
#'   \itemize{
#'     \item modl: model type, i.e., hill
#'     \item fit: fitable?, 1 or 0
#'     \item aic: AIC value
#'     \item tp: model top, <0 means the fit for decreasing direction is preferred
#'     \item ga: ac50
#'     \item gw: Hill coefficient
#'     \item er: scale term
#'    }
#'
#' ## cnst
#'   Fit output from constant model
#'   \itemize{
#'     \item modl: model type, i.e., cnst
#'     \item fit: fittable?, 1 or 0
#'     \item aic: AIC value
#'     \item er: scale term
#'   }
#'
#' @export
#' @seealso [tcpl::tcplObjHill()], [tcpl::tcplObjCnst()], [get_hill_fit_config()]
#'
#' @examples
#'
#' concd <- c(-9, -8, -7, -6, -5, -4)
#' respd <- c(0, 2, 30, 40, 50, 60)
#' maskd <- c(0, 0, 0, 0, 0, 1)
#'
#' # run hill only
#' fit_modls(concd, respd, modls = "hill")
#'
#' # run hill only + increasing direction only
#' fit_modls(concd, respd, modls = "hill", hill_pdir = 1)
#'
#' # run with mask at the highest concentration
#' fit_modls(concd, respd, maskd)
#'
#'
fit_modls <- function(Conc, Resp, Mask = NULL, modls = c("hill", "cnst"), ...) {

  # check the input
  args <- list(...)
  modls <- match.arg(modls, c("hill", "cnst"), several.ok = TRUE)
  args <- .check_modls_args(args, modls)
  outd <- .check_mask_on_concresp(Conc, Resp, Mask)


  new_conc <- outd$Conc
  new_resp <- outd$Resp

  result <- purrr::map(
    modls, ~ do.call(
      fit_modl_in,
      c(list(Conc = new_conc, Resp = new_resp, modl = .x), args)
    )
  ) %>% rlang::set_names(modls)

  # sort the output
  result <- sort_modl_aic(result)
  return(result)
}

#' Fit concentration-response data using one type of models.
#'
#' @inheritParams fit_modls
#' @param modl The model type.
#'
#' @return The respective model output (in a list).
#' @keywords internal
#' @noRd
#'

fit_modl_in <- function(Conc, Resp, modl, ...) {

  l <- list(...)

  # handle model specific configurations
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
      get(eval(func_name)), c(list(Conc = Conc, Resp = Resp), l_p)
    )
  }

  return(result)
}


#' Fit the concentration-response data using the Hill equation.
#'
#' @inheritParams fit_modls
#' @param pdir The preferred direction,
#'   1 is for the increasing direction and -1 is for the decreasing direction.
#'   Default is both values (1, -1) and AIC will be used to select one best fit.
#' @param ... The named input configurations for Hill equation.
#'   This can be used to replace the initial configurations of Hill equation.
#'
#' @seealso [get_hill_fit_config()]
#'
#' @return A list of output parameters from Hill model fit.
#'   If the data is not fittable, the default values for aic, tp, ga, er is NA_real_.
#'   If both directions are used, only the direction with the best fit (by AIC) will be reported.
#'
#' \describe{
#'   \item{modl}{model type, i.e., hill}
#'   \item{fit}{fitable?, 1 or 0}
#'   \item{aic}{AIC value}
#'   \item{tp}{model top, <0 means the fit for decreasing direction is preferred}
#'   \item{ga}{ac50}
#'   \item{gw}{Hill coefficient}
#'   \item{er}{scale term}
#' }
#'
#' @keywords internal
#' @noRd
#'

fit_hill_modl <- function(Conc, Resp, pdir = c(1, -1), ...) {

  args <- list(...)
  if (!all(unique(pdir) %in% c(1, -1))) {rlang::abort("pdir has to be 1 and/or -1")}

  if (length(pdir) > 1) {
    result <- purrr::map(
      pdir, ~ do.call(
        fit_hill_modl_in,
        c(list(Conc = Conc, Resp = Resp, pdir = .x), args)
        )
    )
  } else {
    result <- list(
      do.call(
        fit_hill_modl_in,
        c(list(Conc = Conc, Resp = Resp, pdir = pdir), args)
      )
    )
  }
  result <- select_modl_aic(result)
  return(result)
}

#' Select one out of models based on AIC values.
#'
#'
#' @param lmodls A list of models.
#'
#' @return One model output.
#' @keywords internal
#' @noRd
#'
select_modl_aic <- function(lmodls) {

  if(length(lmodls) > 1) {
    aics <- purrr::map_dbl(lmodls, ~ .x[['aic']])
    if (all(is.na(aics))) {
      result <- lmodls[[1]] # all NAs
    } else {
      result <- lmodls[[which.min(aics)]]
    }
  } else {
    result <- lmodls[[1]]
  }
  return(result)
}

#' Sort a list of models based on AIC values.
#'
#' The model with the lowest AIC value will become the first component of the list.
#'
#'
#' @param lmodls A list of models.
#'
#' @return A sorted list of models.
#' @keywords internal
#' @noRd
#'
sort_modl_aic <- function(lmodls) {

    aics <- purrr::map_dbl(lmodls, ~ .x[['aic']])
    result <- lmodls[order(aics)]

  return(result)
}

#' Fit the concentration-response data using the Hill equation for one direction.
#'
#' @inheritParams fit_hill_modl
#' @param pdir The preferred direction, only one value is allowed.
#'   1 is for the increasing direction and -1 is for the decreasing direction.
#'
#' @inherit fit_hill_modl return
#' @keywords internal
#' @noRd
#'
#'
fit_hill_modl_in <- function(Conc, Resp, pdir = c(1, -1), ...) {

  args <- list(...)

  # default values
  hill <- 0L
  haic <- NA_real_
  hill_tp <- NA_real_
  hill_ga <- NA_real_
  hill_gw <- NA_real_
  hill_er <- NA_real_

  conc <- na.omit(Conc)
  resp <- na.omit(Resp)

  ## handle the cases of all NA ##
  ## special treatment ##
  if (length(resp) == 0) {
    conc <- 0
    resp <- 0
  }

  # adjust the direction
  if (pdir == -1) resp <- resp*-1

  # replace the default config if needed
  config <- get_hill_fit_config(conc, resp)
  config <- .check_hill_args(config, args)

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



#' Get the default configurations for the Hill fit
#'
#' The function gives the default settings by using one set of concentration-response data.
#'
#' @param Conc A vector of log10 concentrations.
#' @param Resp A vector of numeric responses.
#' @param optimf The default optimized function is [tcpl::tcplObjHill()].
#'   but can be changed to ObjHillnorm().
#'
#' @return A list of input configurations.
#'
#' \itemize{
#'   \item theta: initial values of parameters for Hill equation: tp, ga, gw, er
#'   \item f: the object function
#'   \item ui: the bound matrix
#'   \item ci: the bound constraints
#' }
#' @export
#' @seealso [tcpl::tcplObjHill()], [fit_modls()]
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
             logc_min + log10(LBratio), -(logc_max + 0.5), # ga bounds
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



#' Guess the initial scaling factor.
#'
#' @param Resp A vector of resps.
#'
#' @return An initial scaling factor.
#' @keywords internal
#' @noRd

get_est_error <- function(Resp) {
  rmad <- mad(Resp)
  er_est <- if (rmad > 0) log(rmad) else log(1e-32)
  return(er_est)
}



#' Adjust the conc and resp using the mask.
#'
#' @param Conc A vector of log10 concentrations.
#' @param Resp A vector of responses.
#' @param Mask A vector of 1/0 for masking.
#'
#' @return  A list of updated Conc and Resp.
#' @keywords internal
#' @noRd
#'
mask_resp_conc <- function(Conc, Resp, Mask) {

  Resp <- replace(Resp, as.logical(Mask), NA)
  ind <- which(!is.na(Resp))
  Resp <- na.omit(Resp[ind])
  Conc <- na.omit(Conc[ind])

  return(list(Conc = Conc, Resp = Resp))
}


#' Fit responses using a constant model.
#'
#' @param Resp A vector of numeric responses.
#' @param optimf The object function, \code{tcplObjCnst}.
#'
#' @return A list of output parameters from constant model fit.
#'
#' \describe{
#'   \item{modl}{model type, i.e., cnst}
#'   \item{fit}{fitable?, 1 or 0}
#'   \item{aic}{AIC value}
#'   \item{er}{scale term}
#' }
#' @keywords internal
#' @noRd
#' @seealso [tcpl::tcplObjCnst()]
#'
fit_cnst_modl <- function(Resp, optimf = "tcplObjCnst") {

  # default value
  cnst <- 0L
  cnst_er <- NA_real_
  caic <- NA_real_

  resp <- na.omit(Resp)


  # scale term estimate
  er_est <- get_est_error(resp)

  #optimization
  cfit <- optim(er_est,
                fn = get(eval(optimf)),
                method = "Brent",
                lower = er_est - 2,
                upper = er_est + 2,
                control = list(fnscale = -1,
                               reltol = 1e-4,
                               maxit = 500),
                resp = resp)

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
