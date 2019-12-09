
#' @title Default configuration for curvep
#'
#' @description curvep parameters
#'
#' @return a list of parameters with class as curvep_config
#'
#' \itemize{
#'   \item TRSH base(zero-)line threshold
#'   \item RNGE arget range of responses
#'   \item MXDV maximum allowed deviation from monotonicity
#'   \item CARR carryover detection threshold (analysis skipped if set to 0)
#'   \item BSFT for baseline shift issue, min.#points to detect baseline shift (analysis skipped if set to 0)
#'   \item USHP for u-shape curves, min.#points to avoid flattening (analysis skipped if set to 0)
#'   \item TrustHi for equal sets of corrections, trusts those retaining measurements at high concentrations
#'   \item StrictImp prevents extrapolating over concentration-range boundaries; used for POD, ECxx etc.
#'   \item DUMV dummy value for inactive (not suggested to modify)
#'   \item TLOG denominator for calculation wAUC (not suggested to modify)
#' }
#' @export
#' @examples
#'
#' # display all default settings
#' curvep_defaults()
#'
#' # customize settings
#' custom_settings <- curvep_defaults()
#' custom_settings$TRSH <- 30
#' custom_settings
#'
#'
#'
curvep_defaults <- function() {

  defaults <- list(
    TRSH = 15,
    RNGE = -1000000,
    MXDV = 5,
    CARR = 0,
    BSFT = 3,
    USHP = 4,
    TrustHi = TRUE,
    StrictImp = TRUE,
    DUMV = -999,
    TLOG = -24
  )

  class(defaults) <- "curvep_config"
  return(defaults)
}



#' Run rcurvep on a base dataset
#'
#' For a base dataset, only one response is allowed for endpoint-chemical-conc pair.
#' see \code{\link{combi_run_rcurvep} for simulated datasets
#'
#' @param d a dataset with columns: endpoint, chemical, conc, and resp, mask (optional) (see \code{\link{zfishbeh}})
#' @param mask (default = 0, not using mask); if mask column exists, the setting will be ignored
#' use a vector of integers to mask the resps, 1 to mask resp at the highest conc, 2 to mask resp at the second highest conc, and so on.
#' @param config default values are by curvep_defaults()
#' @param keep_sets act_set + (resp_set or fp_set) or all three
#' \itemize{
#'   \item act_set: activity data
#'   \item resp_set: response data
#'   \item fp_set: fingerprint data
#' }
#' @param ... a list of settings; can be used to overwrite the default values in curvep_defaults()
#'
#' @return a rcurvep object; it has two components: result, config
#' @export
#'
#' @examples
#'
#' data(zfishbeh)
#' d <- create_dataset(zfishbeh)
#'
#' # default
#' out <- run_rcurvep(d)
#'
#' # change THR
#' out <- run_rcurvep(d, TRSH = 30)
#'
#' # mask response at highest and second highest concentration
#' out <- run_rcurvep(d, mask = c(1, 2))
#'
#'
run_rcurvep <- function(d, mask = 0, config = curvep_defaults(),  keep_sets = c("act_set", "resp_set", "fp_set"), ...) {

  # not allow NA in the dataset
  d <- na.omit(d)

  d <- .check_dat_base(d)
  mask <- .check_mask_input(mask, d)
  config <- .check_config_name(config = config, ...)
  config <- .check_config_value(config)
  keep_sets <- .check_keep_sets(keep_sets)

  # generate mask column if there is no mask column otherwise pass through the data
  d <- create_resp_mask(d, mask)

  # calculate curvep
  out <- cal_curvep_dataset(d, config = config)

  # clean up the output
  result <- clean_curvep_output(out, config = config)

  # merge the different outputs
  merge_result <- merge_rcurvep_output(d = result, keep_sets = keep_sets)

  out_result <- list(result = merge_result, config = config)
  class(out_result) <- c("rcurvep", class(out_result))

  return(out_result)
}


#' Creat response mask
#'
#' @param d a dataset with columns: endpoint, chemical, conc, and resp, mask (optional) (see \code{\link{zfishbeh}})
#'
#' @param mask 0 for no mask; NULL no change of the data
#' use a vector of integers to mask the resps, 1 to mask resp at the highest conc, 2 to mask resp at the second highest conc, and so on.
#'
#' @return a tibble with added mask column
#' @keywords internal

create_resp_mask <- function(d, mask) {


  if (any(is.null(mask))) { # mask column exists
    result <- d
  } else if (any(mask == 0)) {   # no mask
    result <- d %>% dplyr::mutate(mask = 0)
  } else {

    # generate mask
    result <- d %>%
      dplyr::arrange(.data$endpoint, .data$chemical, dplyr::desc(.data$conc)) %>%
      tidyr::nest(-.data$endpoint, -.data$chemical, .key = "data") %>%
      dplyr::mutate(
        mask = purrr::map(data, function(x, mask) replace(rep(0, nrow(x)), mask, 1), mask = mask)
      ) %>%
      tidyr::unnest()
  }
  return(result)
}

#' Prepare the nested structure of for curvep and calculate (call_curvep)
#'
#' @param d a dataset with columns: endpoint, chemical, conc, and resp, mask
#' @param config curvep_defaults()
#' @return a tibble with a new column output
#' @keywords internal

cal_curvep_dataset <- function(d, config) {
  # prepare the list of data
  d <- d %>%
    dplyr::arrange(.data$endpoint, .data$chemical, .data$conc) %>%
    tidyr::nest(-.data$endpoint, -.data$chemical, .key = "input")

  # use the input and config to call_curvep
  result <- d %>%
    dplyr::mutate(
      output = purrr::map(
        .data$input,
        function(x, config) call_curvep(x$conc, x$resp, x$mask, config), config = config)
    )

  return(result)
}

#' Clean Curvep output to extract the useful information
#'
#' @param d a tibble from cal_curvep_dataset()
#' @param config curvep_defaults()
#'
#' @return a tibble with new columns, out_resp, in_summary, fingerprint, activity
#' @keywords internal

clean_curvep_output <- function(d, config) {
  result <- d %>%
    dplyr::mutate(
      out_resp = purrr::map(.data$output, extract_curvep_outresp),
      in_summary = purrr::map(.data$input, extract_input_summary),
      fingerprint = purrr::map(.data$output, extract_curvep_fingerprint, config = config),
      activity = purrr::map(.data$output, extract_curvep_activity, config = config)
    ) %>%
    dplyr::select(-.data$output)
  return(result)
}


#' Call curvep
#'
#' @param concs a vector of concs (from lowest to highest concentration)
#' @param resps a vector of resps (from lowest to highest concentration)
#' @param masks a vector of masks
#' @param paras curvep_defaults() or a modified list
#' @keywords internal
#' @return (see \code{\link{curvep}}) but remove "Settings"

call_curvep <- function(concs, resps, masks = NULL, paras = curvep_defaults()) {
  result <- do.call(curvep, c(list(concs), list(resps), list(masks), paras))

  # remove the Settings to make list less complicated
  result[['Settings']] <- NULL
  return(result)
}


#' Extract Curvep finterprint
#'
#' @param x the list from curvep output
#' @param config to get the DUMV parameter
#' @keywords internal
#' @return a tibble with three columns: xx, ECxx, Cxx

extract_curvep_fingerprint <- function(x, config) {
  outp <- x
  vals <- outp[names(outp) %in% c('xx', 'ECxx', 'Cxx')]
  result <- vals %>% tibble::as_tibble()

  # make the dummy value as NA
  result[result == config$DUMV] <- NA
  return(result)
}

#' Extract Curvep response output
#'
#' @param x the list from curvep output
#' @keywords internal
#' @return a tibble with two columns: corrected_resp (resp) and correction (corr)

extract_curvep_outresp <- function(x) {
  outp <- x
  vals <- outp[names(outp) %in% c('resp', 'corr')]
  result <- vals %>% tibble::as_tibble()
  result <- result %>%
    dplyr::rename(
      corrected_resp = .data$resp,
      correction = .data$corr
    )
  return(result)
}

#' Extract input (conc, resp) summary data
#'
#' @param x the tibble of the input
#' @keywords internal
#' @return a tibble

extract_input_summary <- function(x) {

  result <- tibble::tibble(
    lowest_conc = head(x$conc, 1),
    highest_conc = tail(x$conc, 1),
    max_resp = max(x$resp),
    min_resp = min(x$resp),
    n_conc = length(unique(x$conc)),
    mean_conc_spacing = mean(diff(unique(x$conc)))

  )
  return(result)
}

#' Extract Curvep activity
#'
#' @param x the list from curvep output
#' @param config for the DUMV parameter
#' @keywords internal
#' @return a tibble with activity
#'
extract_curvep_activity <- function(x, config) {

  outp <- x
  vals <- outp[!names(outp) %in% c('resp', 'corr', 'xx', 'ECxx', 'Cxx','Settings')]
  result <- vals %>% tibble::as_tibble()

  # make the dummy value as NA
  result[result == config$DUMV] <- NA
  return(result)
}

