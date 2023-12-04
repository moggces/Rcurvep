#' Title
#'
#' @param Conc # need to be log10 to be consistent
#' @param Resp
#' @param classSD
#' @param minYrange
#'
#' @return
#' @export
#'
#' @examples
#' concd <- c(-9, -8, -7, -6, -5, -4)
#' respd <- c(0, 2, 30, 40, 50, 60)
#'
#' concd <- c(-11.635385, -9.249779, -7.341294, -5.432809, -4.478566)
#' respd <- c(-9.1669026,  1.2089576, -19.6264789, -56.8283945, -82.5833774)
#'
#'
fit_cc2_modl <- function(Conc, Resp, classSD = 5, minYrange = 20) {

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

  # handle (this has be checked in the input)
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

