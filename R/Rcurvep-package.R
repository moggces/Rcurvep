#' Rcurvep: a R package for concentration response modeling using Curvep
#'
#' The packages provides an easy-to-use R interface for processing concentration response data using Curvep, a response noise filtering algorithm.
#'
#' With a predefined response noise threshold in an experiment, it allows to calculate activity (with confidence interval) based on original or simulated concentration response data.
#' The simulation process differs for three types of datasets:
#' \itemize{
#'   \item binary incidece data (e.g., mortality data)
#'   \item high number of replicates (e.g., alternative animal model data)
#'   \item low number of replicates (e.g., in vitro data)
#' }
#' If the response noise threshold is unknown, the above process can be repeated using a reasonable number of options.
#' The optimal threshold is identified as the one that varience of potency estimation is sufficiently reduced and even stabilized.
#'
#' To learn more about Rcurvep, start with the vignettes:
#' `browseVignettes(package = "Rcurvep")`
#'
#' @import rlang
#' @importFrom magrittr "%>%"
#' @docType package
#' @name Rcurvep
NULL
