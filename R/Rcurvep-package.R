#' Rcurvep: an R package for concentration response modeling using Curvep
#'
#' The packages provides an R interface for processing concentration response data using Curvep, a response noise filtering algorithm.
#'
#' With a predefined baseline noise threshold (or minimum response threshold) in an experiment, it allows to calculate activity (with confidence interval) based on original or simulated concentration response data.
#'
#' If the baseline noise threshold is unknown, the above process can be repeated using a reasonable number of candidates.
#' The optimal threshold is identified as the one that variance of potency estimation is sufficiently reduced and even stabilized.
#'
#' Currently simulated data can be generated from one of the three types of dataset:
#'
#' \itemize{
#'   \item dichotomous binary incidence data (e.g., mortality data from rodent model)
#'   \item continuous data with high number of replicates (e.g., alternative animal model data)
#'   \item continous data with low number of replicates (e.g., in vitro data)
#' }
#'
#' To learn more about Rcurvep, start with the vignettes:
#' `browseVignettes(package = "Rcurvep")`
#'
#' @importFrom magrittr "%>%"
#' @docType package
#' @name Rcurvep
NULL
