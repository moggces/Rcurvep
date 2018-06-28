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
#'
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
#' @seealso \code{\link{curvep}} for available Curvep parameters
#'
#'

extract_curvep_data <- function(c_out, type){

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
    in_concs <- extract_curvep_data(c_out, "concs_in")
    in_resps <- extract_curvep_data(c_out, "resps_in")
    out_resps  <- extract_curvep_data(c_out, "resps_out")

    m <- act %>%
      dplyr::mutate(
        POD = ifelse(is.na(POD), conc_highest, POD),
        EC50 = ifelse(is.na(EC50), conc_highest, EC50)) %>%
      dplyr::group_by(endpoint, chemical, direction, threshold)

    # summarize => confidence interval
    result1 <- m %>%
      dplyr::summarise_at(
        dplyr::vars(dplyr::one_of("POD", "EC50", "Emax", "wAUC", "wAUC_prev")),
        dplyr::funs(med = median(.), ciu = quantile(., probs = 0.975), cil = quantile(., probs = 0.025) )
      ) %>% dplyr::ungroup()

    # summarize => hit confidence
    result2 <- m %>%
      dplyr::summarize(
        hit_confidence = sum(hit)/n()
      ) %>% dplyr::ungroup()

    # summarize => median response in/out
    result3 <- list(in_concs, in_resps, out_resps) %>%
      purrr::reduce(dplyr::inner_join, by = c("repeat_id", "threshold", "endpoint", "chemical", "direction" )) %>%
      tidyr::unnest() %>%
      dplyr::group_by(endpoint, chemical, direction, threshold, concs) %>%
      dplyr::summarize(
        resps = round(median(resps),2),
        resps_in = round(median(resps_in), 2),
      ) %>%
      dplyr::summarize(concs = list(concs), resps = list(resps), resps_in = list(resps_in)) %>%
      dplyr::ungroup()

    result <- list(result1, result2, result3) %>%
      purrr::reduce(dplyr::inner_join, by = c( "threshold", "endpoint", "chemical", "direction" ) )
  }
  return(result)
}
