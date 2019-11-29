
#' Run rcurvep on a base/simulated dataset with a combination of parameters
#'
#' A wrapper function combining \code{\link{create_dataset}}, \code{\link{run_rcurvep}}
#'
#' @param d a dataset such as \code{\link{zfishdev}} and \code{\link{zfishbeh}}
#' @param n_samples NULL (default) or an int to indicate the number of resp per conc to simulate
#' @param vdata NULL (default) or a dbl vector of responses in vehicle control wells; experimental feature
#' @param mask 0 (default) for not using mask
#' use a vector of integers to mask the resps, 1 to mask resp at the highest conc, 2 to mask resp at the second highest conc, and so on.
#' @param keep_sets act_set + (resp_set or fp_set) or all three
#'  \itemize{
#'   \item act_set: activity data
#'   \item resp_set: response data
#'   \item fp_set: fingerprint data
#' }
#' @param ... a list of settings; are used to overwrite the default values in curvep_defaults()
#'
#' @return a list with two components: result (a named list, based on keep_data), config
#' @export
#'
#' @examples
#'
#' data(zfishbeh)
#'
#' out <- combi_run_rcurvep(zfishbeh, n_samples = 2, TRSH = c(5, 10), mask = 1, keep_sets = "act_set")
#'
#'
#'
combi_run_rcurvep <- function(d, n_samples = NULL, vdata = NULL, mask = 0,
                              keep_sets = c("act_set", "resp_set", "fp_set"), ...) {

  paras <- list(...)

  # NA is not allowed
  d <- na.omit(d)

  # check arguments
  d <- .check_dat(d)
  dat_type <- assign_dataset_type(d)
  n_samples <- .check_n_samples(n_samples)
  vdata <- .check_vdata(vdata, dat_type)
  new_config <- .check_config_name2(config = curvep_defaults(), ...)
  keep_sets <- .check_keep_sets(keep_sets)

  # create inputs
  d1 <- create_dataset(d, n_samples = n_samples, vdata = vdata)
  para_in <- create_para_input(paras, n_samples = n_samples, d = d1)

  # run rcurvp using parameters from input
  result <- combi_run_curvep_in(
    d = d1, mask = mask, n_samples = n_samples,
    keep_sets = keep_sets, paras = para_in)

  # unnest the output
  flat_result <- purrr::map(keep_sets, flat_result_tbl, d = result) %>%
    rlang::set_names(keep_sets)

  out_result <- list(result = flat_result, config = new_config)
  class(out_result) <- c("rcurvep", class(out_result))

  return(out_result)

}


#' Merge columns from the result field of an rcurvep object
#'
#' @param rcurvep_obj an object with class:rcurvep
#' @param keep_sets one or multiple of these: act_set, resp_set, fp_set
#'  \itemize{
#'   \item act_set: activity data
#'   \item resp_set: response data
#'   \item fp_set: fingerprint data
#' }
#' @return a list with two components: result (a named list, act_set, resp_set, fp_set), config
#'
#' @keywords internal
#'
merge_rcurvep_output <- function(d, keep_sets) {

  # need to check d
  #rcurvep_obj <- .check_class(rcurvep_obj, "rcurvep", "not a rcurvep object")
  keep_sets <- .check_keep_sets(keep_sets)

  tbl_names <- list(
    act_set = c("in_summary", "activity"),
    resp_set = c("input", "out_resp"),
    fp_set = "fingerprint"
  )

  result <- purrr::map(
    keep_sets, function(x, tbl_names, obj)
      #obj[['result']][c("endpoint", "chemical", tbl_names[[x]])] %>% tidyr::unnest(),
      obj[c("endpoint", "chemical", tbl_names[[x]])] %>% tidyr::unnest(),
    tbl_names = tbl_names,
    obj = d
  ) %>% rlang::set_names(keep_sets)

  #return(list(result = result, config = d$config))
  return(result)
}

#' create parameter input
#'
#' @param paras curvep parameters
#' @param n_samples NULL or an int to indicate the number of resp per conc to simulate
#' @param d a dataset with columns: endpoint, chemical, conc, and resp, mask (optional) (see \code{\link{zfishbeh}})
#'
#' @return a tibble with all combinations of parameters in a data column + sample_id (for a simulated dataset)
#' or a tibble with all combinations of parameters
#' @keywords internal

create_para_input <- function(paras, n_samples, d) {
  if (is.null(n_samples)) {
    result <- expand.grid(paras) %>% tibble::as_tibble()
  } else {
    paras[['sample_id']] <- seq(1, n_samples)
    result <- expand.grid(paras) %>% tibble::as_tibble()
    result <- result %>%
      dplyr::left_join(
        d %>% tidyr::nest(-.data$sample_id, .key = "data"),
        by = "sample_id"
      )
  }
  return(result)
}

#' Run rcurvep using the parameters supplied by purrr::pmap
#'
#' @param d a dataset after \code{\link{create_dataset}}
#' @param mask 0 (default) for not using mask
#' @param n_samples  NULL (default) or an int to indicate the number of resp per conc to simulate
#' @param ... curvep parameters, (and sample_id, data column), created by \code{\link{create_para_input}}
#'
#' @return an rcurvep object
#' @keywords internal
#'
pmap_run_rcurvep <- function(d, mask, n_samples, keep_sets, ...) {

  args <- list(...)

  if (is.null(n_samples)) {
    result <- do.call(run_rcurvep, c(list(d = d, mask = mask, keep_sets = keep_sets), args))
  } else {
    args_f <- args
    args_f[c('data', 'sample_id')] <- NULL
    result <- do.call(run_rcurvep, c(list(d = args$data, mask = mask, keep_sets = keep_sets), args_f))
  }
  return(result)
}

#' Run curvep using created curvep parameters on
#'
#' @param d a dataset
#' @param mask 0 (default) for not using mask
#' @param n_samples NULL (default) or an int to indicate the number of resp per conc to simulate
#' @param paras \code{\link{create_para_input}}
#'
#' @return a tibble (paras) with a new rcurvep_obj column
#' @keywords internal
#'
combi_run_curvep_in <- function(d, mask, n_samples, keep_sets, paras) {
  result <- paras %>%
    dplyr::mutate(
      rcurvep_obj = purrr::pmap(
        ., ~ pmap_run_rcurvep(d = d, mask = mask, n_samples = n_samples, keep_sets = keep_sets, ...)
      )
    )
  suppressWarnings(result <- result %>% dplyr::select(-tidyselect::one_of("data")))
  return(result)
}


#' Map through the rcurvep objects in the tibble and unnest
#'
#' @param d a tibble (paras) with a new rcurvep_obj column \code{\link{combi_run_curvep_in}}
#' @param keep_set only one is allowed
#'
#' @return a tibble
#' @keywords internal
#'
flat_result_tbl <- function(d, keep_set) {
  result <- d %>%
    dplyr::mutate(
      #temp = purrr::map(.data$rcurvep_obj, ~ merge_rcurvep_output(.x, keep_set)$result[[1]])
      temp = purrr::map(.data$rcurvep_obj, ~ .x[['result']][[keep_set]])
    ) %>%
    dplyr::select(-.data$rcurvep_obj) %>%
    tidyr::unnest()

  return(result)

}


