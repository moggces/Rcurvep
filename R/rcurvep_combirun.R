
#' Run Curvep on datasets of concentration-response data with a combination of Curvep parameters
#'
#' It simplifies the steps of [run_rcurvep()] by wrapping the [create_dataset()] in the function.
#'
#' @inheritParams create_dataset
#' @inheritParams run_rcurvep
#' @inherit run_rcurvep return
#' @export
#' @seealso [run_rcurvep()] [summarize_rcurvep_output()]
#'
#' @examples
#'
#' data(zfishbeh)
#'
#' # 2 simulated sample curves +
#' # using two thresholds +
#' # mask the response at the higest concentration
#' # only to output the act_set
#'
#' out <- combi_run_rcurvep(
#'   zfishbeh,
#'   n_samples = 2,
#'   TRSH = c(5, 10),
#'   mask = 1,
#'   keep_sets = "act_set")
#'
#' # create the zfishdev_act dataset
#'\donttest{
#'
#'  data(zfishdev_all)
#'  zfishdev_act <- combi_run_rcurvep(
#'    zfishdev_all, n_samples = 100, keep_sets = c("act_set"),TRSH = seq(5, 95, by = 5),
#'    RNGE = 1000000, CARR = 20, seed = 300
#'  )
#'}
#'
combi_run_rcurvep <- function(d, n_samples = NULL, vdata = NULL, mask = 0,
                              keep_sets = c("act_set", "resp_set", "fp_set"),  ...) {

  paras <- list(...)

  # NA is not allowed
  d <- na.omit(d)

  # check arguments
  d <- .check_dat(d)
  dat_type <- assign_dataset_type(d)
  n_samples <- .check_n_samples(n_samples)
  vdata <- .check_vdata(vdata, dat_type)
  new_config <- .check_config_name2(config = curvep_defaults(), ...)
  keep_sets <- .check_keep_sets(keep_sets, c("act_set", "resp_set", "fp_set"), must_set = "act_set")

  # create inputs
  if (!is.na(new_config$seed)) set.seed(new_config$seed)
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


#' Merge columns from the result field of an rcurvep object.
#'
#' @param rcurvep_obj An object with class:rcurvep.
#' @param keep_sets One or multiple of these sets: act_set, resp_set, fp_set.
#'
#'  \itemize{
#'   \item act_set: activity data
#'   \item resp_set: response data
#'   \item fp_set: fingerprint data
#' }
#' @inherit run_rcurvep return
#' @keywords internal
#' @noRd
#'
merge_rcurvep_output <- function(d, keep_sets) {

  tbl_names <- list(
    act_set = c("in_summary", "activity"),
    resp_set = c("input", "out_resp"),
    fp_set = "fingerprint"
  )

  result <- purrr::map(
    keep_sets, function(x, tbl_names, obj)
      obj[c("endpoint", "chemical", tbl_names[[x]])] %>%
      tidyr::unnest(cols = tbl_names[[x]]),
    tbl_names = tbl_names,
    obj = d
  ) %>% rlang::set_names(keep_sets)

  return(result)
}

#' create Curvep parameter input
#'
#' @param paras Curvep parameters.
#' @inheritParams run_rcurvep
#' @inheritParams create_dataset
#'
#' @return A tibble with all combinations of parameters in a data column + sample_id (for a simulated dataset)
#' or a tibble with all combinations of parameters.
#' @keywords internal
#' @noRd

create_para_input <- function(paras, n_samples, d) {
  if (is.null(n_samples)) {
    result <- expand.grid(paras) %>% tibble::as_tibble()
  } else {
    paras[['sample_id']] <- seq(1, n_samples)
    result <- expand.grid(paras) %>% tibble::as_tibble()
    result <- result %>%
      dplyr::left_join(
        d %>% tidyr::nest(data = -.data$sample_id),
        by = "sample_id"
      )
  }
  return(result)
}

#' Run rcurvep using the parameters supplied by purrr::pmap.
#'
#' @param d Datasets after [create_dataset()].
#' @inheritParams create_dataset
#' @inheritParams run_rcurvep
#' @param ... Curvep parameters, (and sample_id, data column), created by [create_para_input()].
#'
#' @return An rcurvep object.
#' @keywords internal
#' @noRd
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
#' @noRd
#'
combi_run_curvep_in <- function(d, mask, n_samples, keep_sets, paras) {
  result <- paras %>%
    dplyr::mutate(
      rcurvep_obj = furrr::future_pmap(
        ., ~ pmap_run_rcurvep(d = d, mask = mask, n_samples = n_samples, keep_sets = keep_sets, ...)
      )
    )
  suppressWarnings(result <- result %>% dplyr::select(-tidyselect::one_of("data")))
  return(result)
}


#' Map through the rcurvep objects in the tibble and unnest.
#'
#' @param d A tibble (paras) with a new rcurvep_obj column (from [combi_run_curvep_in()]).
#' @param keep_set Only one value is allowed.
#'
#' @return A tibble
#' @keywords internal
#' @noRd
#'
flat_result_tbl <- function(d, keep_set) {
  result <- d %>%
    dplyr::mutate(
      temp = purrr::map(.data$rcurvep_obj, ~ .x[['result']][[keep_set]])
    ) %>%
    dplyr::select(-.data$rcurvep_obj) %>%
    tidyr::unnest("temp")

  return(result)

}


