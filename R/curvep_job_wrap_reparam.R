recreate_curvep_input <- function(input, threshold, direction, paras) {

  rng <- 1000000
  if (direction == -1) { rng <- -1000000 }
  paras['TRSH'] <- threshold
  paras['RNGE'] <- rng
  input[['paras']] <- paras
  return(input)

}



#' Re-paramaterization based on output data frame from run_curvep_job()
#'
#' Given the complex output structure from run_curvep_job(simplify_output = FALSE),
#' the function re-parameterize input for `curvep()`,
#' performs calculations, and generates output.
#'
#' @param dats_out complex output structure from run_curvep_job(simplify_output = FALSE)
#' @param directionality NULL for not changing the directionality, only 1 or -1 is allowed
#' @param threshold NULL for not changing the threshold, only one positive numeric threshold is allowed
#' @param other_paras a list of other Curvep parameters to pass on
#' @param simplify_output (default = FALSE) TRUE to use extract_curvep_data(out, "act")
#'
#' @return see \code{\link{run_curvep_job}}
#' @export
#' @importFrom magrittr "%>%"
#' @importFrom rlang .data
#'
#' @examples
#' data(zfishdev)
#' outd <- run_curvep_job(zfishdev,
#'                       directionality = 1,
#'                       n_sample = 1,
#'                       threshold = 15,
#'                       other_paras = list(CARR = 20, TrustHi = TRUE))
#' outd2 <- reparam_curvep_job(outd, threshold = 25)
#' # more examples are availabie
#' vignette("Rcurvep-intro")
#'
#'
reparam_curvep_job <- function(dats_out, directionality = NULL, threshold = NULL, other_paras = list(), simplify_output = FALSE) {

  #arguments check
  dats_out <- .check_dats_reparam(dats_out)
  directionality <- .check_directionality_reparam(directionality)
  threshold <- .check_threshold_reparam(threshold)
  other_paras <- .check_other_paras(other_paras)


  # change the direction and threshold if needed
  if (!is.null(directionality)) {
    dats_out[, "direction"] <- directionality
  }
  if (!is.null(threshold)) {
    dats_out[, "threshold"] <- threshold
  }

  # special check on the uniqueness
  n_unique <- dats_out %>%
    dplyr::count(threshold, endpoint, chemical, direction, repeat_id) %>%
    dplyr::filter(n == 1) %>%
    nrow()
  if (n_unique != nrow(dats_out)) {
    rlang::abort("only allow unique records (threshold, endpoint, chemical, direction, repeat_id)")
  }

  dats_out <- dats_out %>%
    dplyr::mutate(
      input = purrr::pmap(., function(...) {
        l <- list(...)
        recreate_curvep_input(l$input, l$threshold, l$direction, other_paras)
      })
    )

  #run background curvep
  dats_out <- dats_out %>%
    dplyr::mutate(output = purrr::map(input, run_curvep))

  #get the activity from curvep out
  dats_out <- dats_out %>%
    dplyr::mutate(activity = purrr::map(output, tabulate_curvep_output))

  #simplify the output (especially for BMR finding)
  if (simplify_output) {
    dats_out <- extract_curvep_data(dats_out, "act")
  }

  return(dats_out)
}
