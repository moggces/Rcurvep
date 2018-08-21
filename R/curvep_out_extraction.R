make_comment_as_flat_resp <- function(out_resps, act, comment_pattern) {
  bad_ind <- which(stringr::str_detect(act$Comments, comment_pattern))
  out_resps <- out_resps %>%
    dplyr::mutate(
      resps = dplyr::case_when(
        row_number() %in% bad_ind ~ purrr::map(resps, function(x) {x[x != 0] <- 0; return(x)}),
        TRUE ~ resps
      )
    )
  return(out_resps)
}

make_comment_as_inactive <- function(act, comment_pattern) {
  result <- act %>%
    dplyr::mutate(
      hit = dplyr::case_when(
        stringr::str_detect(Comments, comment_pattern) ~ 0,
        TRUE ~ hit
      ),
      POD = dplyr::case_when(
        stringr::str_detect(Comments, comment_pattern) ~ as.numeric(NA),
        TRUE ~ POD
      ),
      EC50 = dplyr::case_when(
        stringr::str_detect(Comments, comment_pattern) ~ as.numeric(NA),
        TRUE ~ EC50
      ),
      C50 = dplyr::case_when(
        stringr::str_detect(Comments, comment_pattern) ~ as.numeric(NA),
        TRUE ~ C50
      ),
      wConc = dplyr::case_when(
        stringr::str_detect(Comments, comment_pattern) ~ as.numeric(NA),
        TRUE ~ wConc
      ),
      wAUC = dplyr::case_when(
        stringr::str_detect(Comments, comment_pattern) ~ 0,
        TRUE ~ wAUC
      ),
      AUC = dplyr::case_when(
        stringr::str_detect(Comments, comment_pattern) ~ 0,
        TRUE ~ AUC
      ),
      wAUC_prev = dplyr::case_when(
        stringr::str_detect(Comments, comment_pattern) ~ 0,
        TRUE ~ wAUC_prev
      ),
      Emax = dplyr::case_when(
        stringr::str_detect(Comments, comment_pattern) ~ 0,
        TRUE ~ Emax
      ),
      wResp = dplyr::case_when(
        stringr::str_detect(Comments, comment_pattern) ~ 0,
        TRUE ~ wResp
      )
    )
  return(result)
}

make_inactive_potency_as_highest_conc <- function(act) {
  result <- act %>%
    dplyr::mutate(
      POD = ifelse(is.na(POD), conc_highest, POD),
      EC50 = ifelse(is.na(EC50), conc_highest, EC50))
  return(result)
}

add_confidence_interval <- function(m) {
  result <- m %>%
    dplyr::summarise_at(
      dplyr::vars(dplyr::one_of("POD", "EC50", "Emax", "wAUC", "wAUC_prev")),
      dplyr::funs(med = median(.), ciu = quantile(., probs = 0.975), cil = quantile(., probs = 0.025) )
    ) %>% dplyr::ungroup()
  return(result)
}

add_hit_confidence <- function(m) {
  result <- m %>%
    dplyr::summarize(
      hit_confidence = sum(hit)/n()
    ) %>% dplyr::ungroup()
  return(result)
}

calculate_median_resps <- function(in_concs, in_resps, out_resps) {
  result <- list(in_concs, in_resps, out_resps) %>%
    purrr::reduce(dplyr::inner_join, by = c("repeat_id", "threshold", "endpoint", "chemical", "direction" )) %>%
    tidyr::unnest() %>%
    dplyr::group_by(endpoint, chemical, direction, threshold, concs) %>%
    dplyr::summarize(
      resps = round(median(resps),2),
      resps_in = round(median(resps_in), 2),
    ) %>%
    dplyr::summarize(concs = list(concs), resps = list(resps), resps_in = list(resps_in)) %>%
    dplyr::ungroup()
  return(result)
}

#' Unnest Curvep input, output, or activity lists
#'
#' The function unnest the input, output, or activity lists after `run_curvep_job()`
#'
#' @param c_out A tibble after `run_curvep_job()`.
#'
#' @param type A string to indicate which data to extract. Currently seven types are implemented:
#' \itemize{
#'   \item act: all the activity related metrics such as potency and efficacy
#'   \item concs_in: a vector of input concentrations
#'   \item resps_in: a vector of input resps (simulation/original)
#'   \item resps_out: a vector of output (clean) responses
#'   \item paras_in: user manually input parameters
#'   \item paras: all the parameters used in the calculation
#'   \item summary: 1) the hit confidence, 2) the median (med),
#'   95\% confidence interval (ciu, cil) of POD, EC50, Emax, and wAUC, 3) median of resps (input/output)
#' @param modifier A string to match the Curvep Comments column to batch modifiy the activity (active -> inactive)
#' }
#'
#' @return Depending the specified type, a tibble with various columns is returned.
#'
#' @export
#' @importFrom magrittr "%>%"
#' @examples
#' data(zfishdev)
#' x <- zfishdev %>% split(.$endpoint)
#' outd <- run_curvep_job(x[[1]],
#'                       directionality = 1,
#'                       n_sample = 1,
#'                       threshold = 15,
#'                       other_paras = list(CARR = 20, TrustHi = TRUE))
#' # activities
#' extract_curvep_data(outd, "act")
#'
#' # parameters
#' extract_curvep_data(outd, "paras")
#'
#' # summary
#' extract_curvep_data(outd, "summary")
#'
#' # make comments that contain "INVERSE" string as inactive
#'
#' extract_curvep_data(outd, "act", "INVERSE")
#'
#' @seealso \code{\link{curvep}} for available Curvep parameters
#'
#'

extract_curvep_data <- function(c_out, type, modifier = NULL){

  base_ids <- c("repeat_id", "threshold", "endpoint", "chemical", "direction")

  if (type == "act")
  {
    hl <- extract_curvep_data(c_out, "concs_hl")
    result <- c_out %>%
      dplyr::select(dplyr::one_of(base_ids, "activity")) %>%
      tidyr::unnest() %>%
      dplyr::mutate(hit = ifelse(wAUC != 0, 1, 0)) %>%
      dplyr::inner_join(
        hl,  by = c("repeat_id", "threshold", "endpoint", "chemical", "direction")
      )
    if (!is.null(modifier)) {
      result <- make_comment_as_inactive(result, modifier)
    }

  } else if (type == "concs_hl")
  {
    result <- c_out %>%
      dplyr::select(dplyr::one_of(base_ids, "input")) %>%
      dplyr::mutate(
        conc_highest = purrr::map_dbl(input, function(x) max(unlist(x$concs))),
        conc_lowest = purrr::map_dbl(input, function(x) min(unlist(x$concs)))
      ) %>% dplyr::select(-input)
  } else if (type == "concs_in" )
  {
    result <- c_out %>%
      dplyr::select(dplyr::one_of(base_ids, "input")) %>%
      dplyr::mutate(
        concs = purrr::map(input, function(x) unlist(x$concs))
      ) %>% dplyr::select(-input)

  } else if (type == "resps_in" )
  {
    result <- c_out %>%
      dplyr::select(dplyr::one_of(base_ids, "input")) %>%
      dplyr::mutate(
        resps_in = purrr::map(input, function(x) unlist(x$resps))
      ) %>% dplyr::select(-input)

  } else if (type == "resps_out")
  {
    result <- c_out %>%
      dplyr::select(dplyr::one_of(base_ids, "output")) %>%
      dplyr::mutate(resps = purrr::map(output, function(x) x$resp)) %>%
      dplyr::select(-output)

  } else if (type == "paras_in")
  {
    result <- c_out %>%
      dplyr::select(dplyr::one_of(base_ids, "input")) %>%
      dplyr::mutate(temp = purrr::map(paras, function(x) purrr::flatten_df(x$paras))) %>%
      dplyr::select(-input) %>% tidyr::unnest()

  } else if (type == "paras")
  {
    result <- c_out %>%
      dplyr::select(dplyr::one_of(base_ids, "output")) %>%
      dplyr::mutate(temp = purrr::map(output, ~.[['Settings']] %>% purrr::flatten_df())) %>%
      dplyr::select(-output) %>% tidyr::unnest()

  } else if (type == "summary")
  {
    act <- extract_curvep_data(c_out, "act")
    if (!is.null(modifier)) {act <- make_comment_as_inactive(act, modifier)}
    in_concs <- extract_curvep_data(c_out, "concs_in")
    in_resps <- extract_curvep_data(c_out, "resps_in")
    out_resps  <- extract_curvep_data(c_out, "resps_out")
    if (!is.null(modifier)) {
      out_resps <- make_comment_as_flat_resp(out_resps, act, modifier)
    }

    m <- act %>%
      make_inactive_potency_as_highest_conc() %>%
      dplyr::group_by(endpoint, chemical, direction, threshold)

    # summarize => confidence interval ("POD", "EC50", "Emax", "wAUC", "wAUC_prev")
    result1 <- add_confidence_interval(m)

    # summarize => hit confidence
    result2 <- add_hit_confidence(m)

    # summarize => median response in/out
    result3 <- calculate_median_resps(in_concs, in_resps, out_resps)

    result <- list(result1, result2, result3) %>%
      purrr::reduce(dplyr::inner_join, by = c( "threshold", "endpoint", "chemical", "direction" ) )
  }
  return(result)
}
