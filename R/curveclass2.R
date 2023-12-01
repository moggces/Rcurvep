#' Title
#'
#' @param concs # need to be log10 to be consistent
#' @param resps
#' @param sdv
#' @param rnge
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
fit_cc2_modl <- function(concs, resps, sdv = 5, rnge = 10) {

  # create java object and call
  cf <- rJava::.jnew('gov/nih/ncats/ifx/qhts/curvefitting/CurveClassFit3')
  # the input in M
  concs <- 10^concs

  # call the curveclass2
  s <- rJava::.jcall(cf, returnSig = "S", method = "performCurveFitting",
         concs, resps, sdv, rnge)

  # the original output
  #c('ac50', 'curveClass2', 'efficacy', 'hillCoef',
  #  'infActivity', 'logAc50','maskFlag', 'maskNo',
  #  'maxResponse', 'pValue', 'r2', 'zeroActivity')

  # split the output into vector
  # ac50 and maxResponse are not included in the list
  v <- stringr::str_split_fixed(s, pattern = "\n", n = 13)[1:12]

  # handle consistency
  aic <- NA_real_
  fit <- 1
  hillcoef <- as.numeric(v[4])
  if(is.na(hillcoef)) fit <- 0

  # handle cc2 = 4 (no fit)
  masks <- v[7] #null to "0 0 0 0"
  nmasks <- v[8] #NA to 0
  if (v[7] == "null") {
    masks <- stringr::str_c(rep(0, length(concs)), collapse = " ")
    nmasks <- 0
  }

  return(

    list(

      # consistency to tcpl
      aic = aic,
      fit = fit,
      modl = "cc2",

      # original output
      cc2 = as.numeric(v[2]),
      emax = as.numeric(v[3]),
      gw = hillcoef,
      'tp' = as.numeric(v[5]),
      'bt' = as.numeric(v[12]),
      'ga' = as.numeric(v[6]),
      'masks' = masks,
      'nMasks' = as.numeric(nmasks),
      'pvalue' = as.numeric(v[10]),
      'r2' = as.numeric(v[11])

      # added output
      #


    )
  )

}

cc2_2_rank <- function() {

}

#test_fit
#tcplHillVal(inp$conc, win_modl$tp, win_modl$ga, win_modl$gw)
