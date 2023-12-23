#' Fit concentration-response data using Curve Class2 approach
#'
#' Curve Class2 uses 4-parameter Hill model to fit the data.
#' The algorithm assumes the responses are in percentile.
#' Curve Class2 classifies the curves based on fit quality and response magnitude.
#'
#' \describe{
#'   \item{cc2 = 1.1}{2-asymptote curve, pvalue < 0.05, emax > 6\*classSD}
#'   \item{cc2 = 1.2}{2-asymptote curve, pvalue < 0.05, emax <= 6\*classSD & emax > 3\*classSD}
#'   \item{cc2 = 1.3}{2-asymptote curve, pvalue >= 0.05, emax > 6\*classSD}
#'   \item{cc2 = 1.4}{2-asymptote curve, pvalue >= 0.05, emax <= 6\*classSD & emax > 3\*classSD}
#'   \item{cc2 = 2.1}{1-asymptote curve, pvalue < 0.05, emax > 6\*classSD}
#'   \item{cc2 = 2.2}{1-asymptote curve, pvalue < 0.05, emax <= 6\*classSD & emax > 3\*classSD}
#'   \item{cc2 = 2.3}{1-asymptote curve, pvalue >= 0.05, emax > 6\*classSD}
#'   \item{cc2 = 2.4}{1-asymptote curve, pvalue >= 0.05, emax <= 6\*classSD & emax > 3\*classSD}
#'   \item{cc2 = 3}{single point activity, pvalue = NA, emax > 3\*classSD}
#'   \item{cc2 = 4}{inactive, pvalue  >= 0.05, emax <= 3\*classSD}
#'   \item{cc2 = 5}{inconclusive, high bt, further investigation is needed}
#' }
#'
#'
#' @param Conc A vector of log10 concentrations.
#' @param Resp A vector of numeric responses.
#' @param classSD A standard deviation (SD) derived from the responses in the vehicle control.
#'   it is used for classification of the curves. Default = 5%.
#' @param minYrange A minimum response range (max activity - min activity) required to apply curve fitting.
#'   Curve fitting will not be attempted if the response range is less than the cutoff.
#'   Default = 20%.
#' @param ... for additional curve class2 parameters (currently none)
#' @return A list of output parameters from Curve Class2 model fit.
#'   If the data are fit or not fittable (fit = 0), the default value for tp, ga, gw, bt pvalue, masks, nmasks is NA.
#'   For cc2 = 4, it is still possible to have fit parameters.
#'
#'   \itemize{
#'     \item modl: model type, i.e., cc2
#'     \item fit: fittable, 1 (yes) or 0 (no)
#'     \item aic: NA, it is not calculated for this model. The parameter is kept for compatability.
#'     \item cc2: curve class2, default = 4
#'     \item tp: model top, <0 means the fit for decreasing direction is preferred
#'     \item ga: ac50 (log10 scale)
#'     \item gw: Hill coefficient
#'     \item bt: model bottom
#'     \item pvalue: from F-test, for fit quality
#'     \item r2: fitness
#'     \item masks: a string to indicate at which positions of response are masked
#'     \item nmasks: number of masked responses
#'    }
#' @export
#'
#' @references{
#'   \insertRef{PMID:35294762}{Rcurvep}
#' }
#' @seealso  [fit_modls()]
#' @examples
#'
#' fit_cc2_modl(c(-9, -8, -7, -6, -5, -4), c(0, 2, 30, 40, 50, 60))
#'
#'
#'
fit_cc2_modl <- function(Conc, Resp, classSD = 5, minYrange = 20, ...) {

  args <- list(...)
  args <- .check_modl_args_exist(list(NULL), args)
  # the input in M
  Conc <- 10^Conc
  conc <- na.omit(Conc)
  resp <- na.omit(Resp)

  # default values
  cc2_fit <- 0L
  cc2_cc2 <- 4
  cc2_aic <- NA_real_
  cc2_tp <- NA_real_
  cc2_ga <- NA_real_
  cc2_gw <- NA_real_
  cc2_bt <- NA_real_
  cc2_pvalue <- NA_real_
  cc2_r2 <- NA_real_
  cc2_masks <- NA_character_
  cc2_nmasks <- NA_integer_

  ## handle the cases of all NA ##
  ## special treatment ##
  if (length(resp) <= 1) { # it has to be an array to fullfill java signature
    return(
      list(
        # consistency to tcpl
        aic = cc2_aic,
        fit = cc2_fit,
        modl = "cc2",

        # original output
        cc2 = cc2_cc2,
        gw = cc2_gw,
        tp = cc2_tp,
        bt = cc2_bt,
        ga = cc2_ga,
        masks = cc2_masks,
        nMasks = cc2_nmasks,
        pvalue = cc2_pvalue,
        r2 = cc2_r2
      )
    )
  }

  # handle different number of conc and resp
  if(length(conc) != length(resp)) {
    rlang::abort("The length of conc and resp is not the same after removing NA.")
  }

  # create java object and call
  cf <- rJava::.jnew('gov/nih/ncats/ifx/qhts/curvefitting/CurveClassFit3')

  # call the curveclass2
  s <- rJava::.jcall(cf, returnSig = "S", method = "performCurveFitting",
                     conc, resp, classSD, minYrange)

  # the original output
  #c('ac50', 'curveClass2', 'efficacy', 'hillCoef',
  #  'infActivity', 'logAc50','maskFlag', 'maskNo',
  #  'maxResponse', 'pValue', 'r2', 'zeroActivity')

  # split the output into vector
  # ac50 and maxResponse and efficacy are not included in the list
  v <- stringr::str_split_fixed(s, pattern = "\n", n = 13)[1:12]

  # get hill coef first
  suppressWarnings( #null to NA
    hillcoef <- as.numeric(v[4])
  )
  # NULL cannot be returned in the ifelse
  if(!is.na(hillcoef)) cc2_masks <- v[7]


  return(
    suppressWarnings( #null to NA
      list(

        # consistency to tcpl
        aic = cc2_aic,
        fit = ifelse(!is.na(hillcoef), 1L, cc2_fit),
        modl = "cc2",

        # original output
        cc2 = as.numeric(v[2]),
        gw = hillcoef,
        tp = as.numeric(v[5]),
        bt = as.numeric(v[12]),
        ga = as.numeric(v[6]),
        masks = cc2_masks,
        nMasks = ifelse(!is.na(hillcoef), as.numeric(v[8]), cc2_nmasks),
        pvalue = as.numeric(v[10]),
        r2 = as.numeric(v[11])

        # added output
        #

      )
    )
  )

}

cc2_2_rank <- function() {

}

