
#' Run rcurvep on a dataset with a combination of parameters
#'
#' @param d datasets such as \code{\link{zfishdev}} and \code{\link{zfishbeh}}
#' @param n_samples NULL (default) or an int to indicate the number of resp per conc to simulate
#' @param vdata NULL (default) or a dbl vector of responses in vehicle control wells
#' @param mask NULL to not use mask to mask resps (default = NULL);
#' use a vector of integers to mask the resps, 1 to mask resp at the highest conc, 2 to mask resp at the second highest conc, and so on.
#' @param keep_data allowed one or multiple values: act_set, resp_set, fp_set
#' @param ... a list of settings; are used to overwrite the default values in curvep_defaults()
#'
#' @return a list with two components: result (a named list, based on keep_data), config
#' @export
#'
#' @examples
#'
#'
combi_run_rcurvep <- function(d, n_samples = NULL, vdata = NULL, mask = NULL,
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
  result <- para_in %>%
    dplyr::mutate(
      rcurvep_obj = purrr::pmap(
        ., ~ pmap_run_rcurvep(d = d1, mask = mask, n_samples = n_samples, ...)
        )
    )
  suppressWarnings(result <- result %>% dplyr::select(-tidyselect::one_of("data")))

  # unnest the output
  flat_result <- purrr::map(
    keep_sets, function(x, result) flat_result_tbl(result, x), result = result
    ) %>% rlang::set_names(keep_sets)

  return(list(result = flat_result, config = new_config))

}


#' Merge columns from the result field of an rcurvep object
#'
#' @param rcurvep_obj an object with class:rcurvep
#' @param keep_sets allowed one or multiple values: act_set, resp_set, fp_set
#'
#' @return a list with two components: result (a named list, based on keep_sets), config
#' @export
#'
merge_rcurvep_output <- function(rcurvep_obj, keep_sets = c("act_set", "resp_set", "fp_set")) {

  rcurvep_obj <- .check_class(rcurvep_obj, "rcurvep", "not a rcurvep object")
  keep_sets <- .check_keep_sets(keep_sets)

  tbl_names <- list(
    act_set = c("in_summary", "activity"),
    resp_set = c("input", "out_resp"),
    fp_set = "fingerprint"
  )

  result <- purrr::map(
    keep_sets, function(x, tbl_names, obj)
      obj[['result']][c("endpoint", "chemical", tbl_names[[x]])] %>% tidyr::unnest(),
    tbl_names = tbl_names,
    obj = rcurvep_obj
  ) %>% rlang::set_names(keep_sets)

  return(list(result = result, config = rcurvep_obj$config))
}

#' create parameter input
#'
#' @param paras curvep parameters
#' @param n_samples
#' @param d
#'
#' @return a tibble with all combinations of parameters
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
#' @param d
#' @param mask
#' @param n_samples
#' @param ... curvep parameters, (and sample_id, data column)
#'
#' @return an rcurvep object
#' @keywords internal
#'
pmap_run_rcurvep <- function(d, mask, n_samples, ...) {

  args <- list(...)

  if (is.null(n_samples)) {
    result <- do.call(run_rcurvep, c(list(d = d, mask = mask), args))
  } else {
    args_f <- args
    args_f[c('data', 'sample_id')] <- NULL
    result <- do.call(run_rcurvep, c(list(d = args$data, mask = mask), args_f))
  }
  return(result)
}


#' Map through the rcurvep objects in the tibble and unnest
#'
#' @param d
#' @param keep_set only one is allowed
#'
#' @return a tibble
#' @keywords internal
#'
flat_result_tbl <- function(d, keep_set) {
  result <- d %>%
    dplyr::mutate(
      temp = purrr::map(.data$rcurvep_obj, ~ merge_rcurvep_output(.x, keep_set)$result[[1]])
    ) %>%
    dplyr::select(-.data$rcurvep_obj) %>%
    tidyr::unnest()

  return(result)

}


