#' Default configuration for curvep
#'
#' curvep_defaults() return a list of parameters with class as "curvep_config"
#'
#' @return a list of parameters with class as "curvep_config"
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
#'
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
#' @export
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



#' Run rcurvep on a dataset
#'
#'
#' @param d a dataset with columns: endpoint, chemical, conc, and resp, mask (optional) (see \code{\link{zfishbeh})
#' @param mask NULL to not use mask to mask resps (default = NULL);
#' use a vector of integers to mask the resps, 1 to mask resp at the highest conc, 2 to mask resp at the second highest conc, and so on.
#' @param config default values are by curvep_defaults()
#' @param ... a list of settings; can be used to overwrite the default values in curvep_defaults()
#'
#' @return a rcurvep object; it has two components: result, config.
#' Result has columns: endpoint, chemical, input, out_resp, in_summary, fingerprint, activity
#'
#' @examples
#'
#'
#'
#' @importFrom rlang .data
#' @export
#'
#'
run_rcurvep <- function(d, mask = NULL, config = curvep_defaults(), ...) {

  # not allow NA in the dataset
  d <- na.omit(d)

  # check data
  d <- .check_dat_base(d)
  mask <- .check_mask_input(mask, d)
  config <- .check_config_name(config = config, ...)
  config <- .check_config_value(config)

  # generate mask column
  d <- create_mask(d, mask)

  # prepare the list of data
  d <- d %>%
    dplyr::arrange(.data$endpoint, .data$chemical, .data$conc) %>%
    tidyr::nest(-.data$endpoint, -.data$chemical, .key = "input")

  # generate output and activity column
  result <- d %>%
    dplyr::mutate(
      output = purrr::map(.data$input, function(x, config) call_curvep(x$conc, x$resp, x$mask, config), config = config),
      out_resp = purrr::map(.data$output, extract_curvep_outresp),
      in_summary = purrr::map(.data$input, extract_input_summary),
      fingerprint = purrr::map(.data$output, extract_curvep_fingerprint, config = config),
      activity = purrr::map(.data$output, extract_curvep_activity, config = config)
    ) %>%
    dplyr::select(-.data$output)

  out_result <- list(result = result, config = config)
  class(out_result) <- 'rcurvep'

  return(out_result)
}


#' Creat mask
#'
#' @return a tibble/data.frame with added mask column
#' @keywords internal
#' @importFrom rlang .data
#'
#'
create_mask <- function(d, mask) {

  result <- d
  # if NULL use all 0 to avoid
  if (is.null(mask)) {
    if (!rlang::has_name(d, "mask")) {
      result <- d %>%
        dplyr::mutate(
          mask = 0
        )
    }
  } else {
    if (rlang::has_name(d, "mask")) d$mask <- NULL
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


#' Call curvep
#'
#' @param concs a vector of concs (from lowest to highest concentration)
#' @param resps a vector of resps (from lowest to highest concentration)
#' @param masks a vector of masks
#' @param paras
#' @keywords internal
#' @return (see \code{\link{curvep}) but remove "Settings"

call_curvep <- function(concs, resps, masks = NULL, paras = curvep_defaults()) {
  result <- do.call(curvep, c(list(concs), list(resps), list(masks), paras))

  # remove the Settings to make list less complicated
  result[['Settings']] <- NULL
  return(result)
}


#' extract_curvep_fingerprint
#'
#' @param x the list from curvep output
#' @param config
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

#' extract_curvep_outresp
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

#' extract_input_summary
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
    n_conc = length(x$conc),
    mean_conc_spacing = mean(diff(x$conc))

  )
  return(result)
}

#' extract_curvep_activity
#'
#' @param x the list from curvep output
#' @param config
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

