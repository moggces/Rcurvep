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

add_confidence_interval <- function(m, col_names = c("POD", "EC50", "Emax", "wAUC", "wAUC_prev"), level = 0.95) {
  upper_bound <- 1 - (1 - level)/2
  lower_bound <- (1 - level)/2
  result <- m %>%
    dplyr::summarise_at(
      dplyr::vars(dplyr::one_of(col_names)),
      dplyr::funs(med = median(.), ciu = quantile(., probs = upper_bound), cil = quantile(., probs = lower_bound) )
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

  base_ids <- c("threshold", "endpoint", "chemical", "direction")
  result <- list(in_concs, in_resps, out_resps) %>%
    purrr::reduce(dplyr::inner_join, by = c("repeat_id", base_ids )) %>%
    tidyr::unnest() %>%
    dplyr::group_by_at(dplyr::vars(dplyr::one_of(c(base_ids, "concs")))) %>%
    dplyr::summarize(
      resps = round(median(resps),2),
      resps_in = round(median(resps_in), 2),
    ) %>%
    dplyr::summarize(concs = list(concs), resps = list(resps), resps_in = list(resps_in)) %>%
    dplyr::ungroup()
  return(result)
}

#' Unnest Curvep input, output, or activity lists.
#'
#' The function unnest the input, output, or activity lists after `run_curvep_job()`.
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
#' }
#' @param modifier A string to match the Curvep Comments column to batch modifiy the activity,
#' including hit, potency related values, and adjusted responses (active -> inactive)
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

  # test the class
  if (!is(c_out, "rcurvep_out_nested")) {
    warning("input is not the nested data frame from the run_curvep_job()")
  }

  base_ids <- c("threshold", "endpoint", "chemical", "direction")

  if (type == "act")
  {
    hl <- extract_curvep_data(c_out, "concs_hl")
    result <- c_out %>%
      dplyr::select(dplyr::one_of(base_ids, "repeat_id", "activity")) %>%
      tidyr::unnest() %>%
      dplyr::mutate(hit = ifelse(wAUC != 0, 1, 0)) %>%
      dplyr::inner_join(
        hl,  by = c(base_ids, "repeat_id")
      )
    if (!is.null(modifier)) {
      result <- make_comment_as_inactive(result, modifier)
    }

  } else if (type == "concs_hl")
  {
    result <- c_out %>%
      dplyr::select(dplyr::one_of(base_ids, "repeat_id", "input")) %>%
      dplyr::mutate(
        conc_highest = purrr::map_dbl(input, function(x) max(unlist(x$concs))),
        conc_lowest = purrr::map_dbl(input, function(x) min(unlist(x$concs)))
      ) %>% dplyr::select(-input)
  } else if (type == "concs_in" )
  {
    result <- c_out %>%
      dplyr::select(dplyr::one_of(base_ids, "repeat_id", "input")) %>%
      dplyr::mutate(
        concs = purrr::map(input, function(x) unlist(x$concs))
      ) %>% dplyr::select(-input)

  } else if (type == "resps_in" )
  {
    result <- c_out %>%
      dplyr::select(dplyr::one_of(base_ids, "repeat_id", "input")) %>%
      dplyr::mutate(
        resps_in = purrr::map(input, function(x) unlist(x$resps))
      ) %>% dplyr::select(-input)

  } else if (type == "resps_out")
  {
    result <- c_out %>%
      dplyr::select(dplyr::one_of(base_ids, "repeat_id", "output")) %>%
      dplyr::mutate(resps = purrr::map(output, function(x) x$resp)) %>%
      dplyr::select(-output)
    if (!is.null(modifier)) {
      result <- make_comment_as_flat_resp(result, act, modifier)
    }

  } else if (type == "paras_in")
  {
    result <- c_out %>%
      dplyr::select(dplyr::one_of(base_ids, "repeat_id", "input")) %>%
      dplyr::mutate(temp = purrr::map(paras, function(x) purrr::flatten_df(x$paras))) %>%
      dplyr::select(-input) %>% tidyr::unnest()

  } else if (type == "paras")
  {
    result <- c_out %>%
      dplyr::select(dplyr::one_of(base_ids, "repeat_id", "output")) %>%
      dplyr::mutate(temp = purrr::map(output, ~.[['Settings']] %>% purrr::flatten_df())) %>%
      dplyr::select(-output) %>% tidyr::unnest()

  } else if (type == "summary")
  {
    act <- extract_curvep_data(c_out, "act")
    in_concs <- extract_curvep_data(c_out, "concs_in")
    in_resps <- extract_curvep_data(c_out, "resps_in")
    out_resps  <- extract_curvep_data(c_out, "resps_out")
    sum_input <- list(act, in_concs, in_resps, out_resps) %>%
      purrr::reduce(dplyr::inner_join, by = c(base_ids, "repeat_id"))

    result <- summarize_curvep_output(sum_input, col_names = c("POD", "EC50", "Emax", "wAUC", "wAUC_prev"), modifier = modifier, conf_level = 0.95)
  }

  # append the class information
  if (type != "summary") {
    class(result) <- c("rcurvep_out", class(result))
  }
  return(result)
}

#' Summarize the activity output after curvep based on bootstrap samples.
#' The summarized values include hit_confidence, confidence interval for the activity values,
#' and median responses per concentration.
#'
#' @param a_out a tibble from the extract_curvep_data(type != "summary") or run_curvep_job(simplify_output = TRUE)
#' @param col_names default values include POD, EC50, Emax, wAUC, wAUC_prev
#' @param modifier A string to match the Curvep Comments column to batch modifiy the activity,
#' including hit, potency related values, and adjusted responses (active -> inactive)
#' @param conf_level default value = 0.95 (95\%)
#'
#' @return a tibble with summary data: 1) the hit confidence, 2) the median (med),
#'   95\% confidence interval (ciu, cil) of activity values, 3) median of responses of concentrations (input/output)
#' 3) is only calculated when input tables have resps, concs, and resps_in columns
#' @export
#'
#' @examples
#' data(zfishdev)
#' x <- zfishdev %>% split(.$endpoint)
#' acts <- run_curvep_job(x[[1]],
#'                       directionality = 1,
#'                       n_sample = 10,
#'                       threshold = 15,
#'                       other_paras = list(CARR = 20, TrustHi = TRUE), simplify_output = TRUE)
#' # activities
#' acts_sum <- summarize_curvep_output(acts, conf_level = 0.9)
#'
summarize_curvep_output <- function(a_out, col_names = c("POD", "EC50", "Emax", "wAUC", "wAUC_prev"), modifier = NULL, conf_level = 0.95) {


  # test the class
  if (!is(a_out, "rcurvep_out")) {
    warning("input is not the activity file from extract_curvep_data() ")
  }

  base_ids <- c("threshold", "endpoint", "chemical", "direction")

  #adjust the potency based on comments
  if (!is.null(modifier)) {a_out <- make_comment_as_inactive(a_out, modifier)}

  #adjust the potency of inactives
  m <- a_out %>%
    make_inactive_potency_as_highest_conc() %>%
    dplyr::group_by_at(dplyr::vars(dplyr::one_of(c(base_ids))))

  # summarize => confidence interval ("POD", "EC50", "Emax", "wAUC", "wAUC_prev")
  result1 <- add_confidence_interval(m, col_names = col_names, level = conf_level)

  # summarize => hit confidence
  result2 <- add_hit_confidence(m)

  # summarize => median response in/out
  if (sum(rlang::has_name(a_out, c("resps", "concs", "resps_in")))  == 3) {


    in_concs <- a_out %>% dplyr::select(dplyr::one_of(base_ids, "repeat_id", "concs"))
    out_resps <- a_out %>% dplyr::select(dplyr::one_of(base_ids, "repeat_id", "resps"))
    in_resps <- a_out %>% dplyr::select(dplyr::one_of(base_ids, "repeat_id", "resps_in"))
    if (!is.null(modifier)) {
      out_resps <- make_comment_as_flat_resp(out_resps, a_out, modifier)
    }
    result3 <- calculate_median_resps(in_concs, in_resps, out_resps)
    result <- list(result1, result2, result3) %>%
      purrr::reduce(dplyr::inner_join, by = c(base_ids) )
  } else {
    result <- list(result1, result2) %>%
      purrr::reduce(dplyr::inner_join, by = c(base_ids) )
  }

  return(result)

}
