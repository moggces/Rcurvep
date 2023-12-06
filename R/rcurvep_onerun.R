
#' Default parameters of Curvep
#'
#' @seealso [curvep()]
#'
#' @return A list of parameters with class as curvep_config.
#'
#' \itemize{
#'   \item TRSH: (default = 15) base(zero-)line threshold
#'   \item RNGE: (default = -1000000, decreasing) target range of responses
#'   \item MXDV: (default = 5) maximum allowed deviation from monotonicity
#'   \item CARR: (default = 0) carryover detection threshold (analysis skipped if set to 0)
#'   \item BSFT: (default = 3) for baseline shift issue, min.#points to detect baseline shift (analysis skipped if set to 0)
#'   \item USHP: (default = 4) for u-shape curves, min.#points to avoid flattening (analysis skipped if set to 0)
#'   \item TrustHi: (default = TRUE)for equal sets of corrections, trusts those retaining measurements at high concentrations
#'   \item StrictImp: (default = TRUE) prevents extrapolating over concentration-range boundaries; used for POD, ECxx etc.
#'   \item DUMV: (default = -999) dummy value for inactive (not suggested to modify)
#'   \item TLOG: (default = -24) denominator for calculation wAUC (not suggested to modify)
#'   \item seed: (default = NA) can be set when bootstrapping samples
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
    TLOG = -24,
    seed = NA_integer_
  )

  class(defaults) <- "curvep_config"
  return(defaults)
}



#' Run Curvep on datasets of concentration-response data
#'
#'
#' The concentration-response relationship per endpoint and chemical has to be 1-to-1.
#' If not, use [create_dataset()] for pre-processing or
#' use [combi_run_rcurvep()], which has both pre-processing and more flexible parameter controls.
#'
#'
#' @param d Datasets with columns: endpoint, chemical, conc, and resp, mask (optional)
#'   Example datasets as [zfishbeh].
#'   It is required that the baseline of responses in the resp column to be 0.
#' @param mask Default = 0, for no mask (values in the mask column all 0).
#'   Use a vector of integers to mask the responses:
#'   1 to mask the response at the highest concentration;
#'   2 to mask the response at the second highest concentration, and so on.
#'   If mask column exists, the setting will be ignored.
#'
#' @param config Default configurations set by [curvep_defaults()].
#' @param keep_sets The types of output to be reported.
#'   Allowed values: act_set, resp_set, fp_set. Multiple values are allowed.
#'   act_set is the must.
#' \itemize{
#'   \item act_set: activity data
#'   \item resp_set: response data
#'   \item fp_set: fingerprint data
#' }
#' @param ... Curvep settings.
#'   See [curvep_defaults()] for allowed parameters.
#'   These can be used to overwrite the default values.
#'
#' @return An rcurvep object. It has two components: result, config
#'   The result component is also a list of output sets depending on the parameter, *keep_sets*.
#'   The config component is a *curvep_config* object.\cr
#'
#'   Often used columns in the *act_set*: AUC (area under the curve), wAUC (weighted AUC),
#'   POD (point-of-departure), EC50 (Half maximal effective concentration),
#'   nCorrected (number of corrected points).
#'
#'
#'
#' @seealso [create_dataset()], [combi_run_rcurvep()], [curvep_defaults()].
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
#' # change TRSH
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
  keep_sets <- .check_keep_sets(keep_sets, c("act_set", "resp_set", "fp_set"), must_set = "act_set")

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


#' Creat a mask column on the dataset.
#'
#' @inheritParams run_rcurvep
#' @return A tibble with a mask column.
#' @keywords internal
#' @noRd

create_resp_mask <- function(d, mask) {


  if (any(is.null(mask))) { # mask column exists
    result <- d
  } else if (any(mask == 0)) {   # no mask
    result <- d %>% dplyr::mutate(mask = 0)
  } else {

    # generate mask
    result <- d %>%
      dplyr::arrange(endpoint, chemical, dplyr::desc(conc)) %>%
      tidyr::nest(data = -c(endpoint, chemical)) %>%
      dplyr::mutate(
        mask = purrr::map(data, function(x, mask) replace(rep(0, nrow(x)), mask, 1), mask = mask)
      ) %>%
      tidyr::unnest(cols = c("data", "mask"))
  }
  return(result)
}

#' Prepare the nested structure of for curvep and calculate (call_curvep).
#'
#' @param d Datasets with columns: endpoint, chemical, conc, resp, and mask.
#' @param config Updated curvep parameters.
#' @return A tibble with a new column, output.
#' @keywords internal
#' @noRd

cal_curvep_dataset <- function(d, config) {
  # prepare the list of data
  d <- d %>%
    dplyr::arrange(endpoint, chemical, conc) %>%
    tidyr::nest(input = -c(endpoint, chemical))

  # use the input and config to call_curvep
  result <- d %>%
    dplyr::mutate(
      output = purrr::map(
        input,
        function(x, config) call_curvep(x$conc, x$resp, x$mask, config), config = config)
    )

  return(result)
}

#' Clean Curvep output to extract the useful information.
#'
#' @param d A tibble from [cal_curvep_dataset()].
#' @inheritParams cal_curvep_dataset
#' @return A tibble with new columns, out_resp, in_summary, fingerprint, activity.
#' @keywords internal
#' @noRd

clean_curvep_output <- function(d, config) {
  result <- d %>%
    dplyr::mutate(
      out_resp = purrr::map(output, extract_curvep_outresp),
      in_summary = purrr::map(input, extract_input_summary),
      fingerprint = purrr::map(output, extract_curvep_fingerprint, config = config),
      activity = purrr::map(output, extract_curvep_activity, config = config)
    ) %>%
    dplyr::select(-output)
  return(result)
}


#' Call curvep function.
#'
#' @seealso [curvep()]
#'
#' @param concs A vector of concentrations (from lowest to highest).
#' @param resps A vector of responses (matching up the concs).
#' @param masks A vector of masks.
#' @param paras Updated curvep parameters.
#' @return The output from [curvep()] but remove "Settings".
#' @keywords internal
#' @noRd

call_curvep <- function(concs, resps, masks, paras) {
  result <- do.call(curvep, c(list(concs), list(resps), list(masks), paras))

  # remove the Settings to make list less complicated
  result[['Settings']] <- NULL
  return(result)
}


#' Extract Curvep finterprint.
#'
#' @param x The list from curvep output.
#' @param config The curvep configuration information (to get the DUMV setting).
#' @return A tibble with three columns: xx, ECxx, Cxx
#' @keywords internal
#' @noRd

extract_curvep_fingerprint <- function(x, config) {
  outp <- x
  vals <- outp[names(outp) %in% c('xx', 'ECxx', 'Cxx')]
  result <- vals %>% tibble::as_tibble()

  # make the dummy value as NA
  result[result == config$DUMV] <- NA
  return(result)
}

#' Extract Curvep response output.
#'
#' @inheritParams extract_curvep_fingerprint
#' @return A tibble with two columns: corrected_resp (resp) and correction (corr)
#' @keywords internal
#' @noRd
#'

extract_curvep_outresp <- function(x) {
  outp <- x
  vals <- outp[names(outp) %in% c('resp', 'corr')]
  result <- vals %>% tibble::as_tibble()
  result <- result %>%
    dplyr::rename(
      corrected_resp = resp,
      correction = corr
    )
  return(result)
}

#' Extract input (conc, resp) summary data.
#'
#' @param x The tibble of the input.
#' @keywords internal
#' @return A tibble.
#' @noRd

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

#' Extract Curvep activity.
#'
#' @inheritParams extract_curvep_fingerprint
#' @keywords internal
#' @return A tibble with activity related data.
#' @noRd
#'
extract_curvep_activity <- function(x, config) {

  outp <- x
  vals <- outp[!names(outp) %in% c('resp', 'corr', 'xx', 'ECxx', 'Cxx','Settings')]
  result <- vals %>% tibble::as_tibble()

  # make the dummy value as NA
  result[result == config$DUMV] <- NA
  return(result)
}

