#' Clean and summarize the output of rcurvep object
#'
#' @details
#' The function can perform the following tasks:
#' \enumerate{
#'   \item add an column, hit, in the *act_set*
#'   \item unhit (make result as inactive) if the Comments column contains a certain string
#'   \item summarize the results
#'}
#' The curve is considered as "hit" if its responses are monotonic after processing by Curvep.
#' However, often, if the curve is "INVERSE" (yet monotonic) is not considered as an active curve.
#' By using the information in the Comments column, we can "unhit" these cases.
#'
#' When (clean_only = FALSE, default), a tibble, act_summary is generated with confidence intervals of the activity metrics.
#' The quantile approach is used to calculate the confidence interval.
#' For potency activity metrics, if value is NA, highest tested concentration is used in the summary.
#' For other activity metrics, if value is NA, 0 is used in the summary.
#'
#'
#' @param d The rcurvep object from [combi_run_rcurvep()] and [run_rcurvep()].
#' @param inactivate A character string, default = NULL,
#'   to make the curve with this string in the Comments column as inactive.
#'   or a vector of index for the rows in the act_set that needs to be inactive
#' @param ci_level Default = 0.95 (95 percent of confidence interval).
#' @param clean_only Default = FALSE, only the 1st, 2nd task will be performed (see Details).
#'
#' @return A list of named components: result and config (and act_summary).
#'   The result and config are the copy of the input d (but with modifications if *inactivate* is not NULL).
#'   If (clean_only = FALSE), an *act_summary* is added.
#'
#'   Suffix meaning in column names in *act_summary*: med (median), cil (lower end confidence interval),
#'   ciu (higher end confidence interval)
#'   Often used columns in *act_summary*: n_curves (number of curves used in summary), hit_confidence (fraction of active in n_curves)
#'
#' @export
#' @seealso [combi_run_rcurvep()], [run_rcurvep()]
#' @examples
#'
#' data(zfishbeh)
#'
#' # original datasets
#' out <- combi_run_rcurvep(zfishbeh, n_samples = NULL, TRSH = c(5, 10))
#' out_res <- summarize_rcurvep_output(out)
#'
#'\dontrun{
#' # unhit when comment has "INVERSE"
#' out <- summarize_rcurvep_output(out, inactivate = "INVERSE")
#'
#' # unhit for certain rows in act_set
#' out <- summarize_rcurvep_output(out, inactivate = c(2,3))
#'
#' # simulated datasets
#' out <- combi_run_rcurvep(zfishbeh, n_samples = 3, TRSH = c(5, 10))
#' out_res <- summarize_rcurvep_output(out)
#'}
#'
summarize_rcurvep_output <- function(d, inactivate = NULL, ci_level = 0.95, clean_only = FALSE) {

  inactivate <- .check_inactivate(inactivate)
  d <- .check_combirun_out(d)
  lsets <- d$result
  lsets <- .check_result_sets(lsets)

  # get the common column names (e.g., TRSH)
  base_cols <- get_base_cols(lsets)

  # add the hit column in the act_set
  lsets$act_set <- add_hit_actset(lsets$act_set)
  lsets$act_set <- apply_comment_to_unhit_actset(lsets$act_set, inactivate = inactivate)
  lsets$act_set <- adjust_rcurvep_comment(lsets$act_set, inactivate = inactivate)

  # nest the sets into columns
  lsets_n <- get_nested_joined_sets(lsets, base_cols)

  # use the Comments column to unhit
  lsets_n <- apply_comment_to_unhit(lsets_n, inactivate = inactivate)

  # unnest the clean sets
  lsets_clean <- purrr::map(
    names(lsets), unnest_joined_sets, nested = lsets_n, base_cols = base_cols) %>%
    rlang::set_names(names(lsets))

  # make the result as list
  # two components: result, config
  result <- list(
    result = lsets_clean,
    config = d$config
  )

  # summarize the act_set
  # three components: result, act_summary, config
  if(!clean_only) {
    # make the activity column = NA as the highest conc
    lsets_n <- make_act_na_highconc(lsets_n)
    lsets_n <- summarize_actsets(lsets_n, ci_level = ci_level)
    result[['act_summary']] <- unnest_joined_sets(lsets_n, base_cols, "act_summary")

  }

  class(result) <- c("rcurvep", class(result))


  return(result)
}


#' Get common column names between sets.
#'
#' Depending on different settings (e.g., TRSH), the common column names change.
#'
#' @param lsets The result list from the [combi_run_rcurvep()] or [run_rcurvep()].
#' @param config curvep_defaults()
#' @param remove_sample_id default = TRUE
#'
#' @return A vector of common column names in sets.
#' @keywords internal
#' @noRd

get_base_cols <- function(lsets, config = curvep_defaults(), remove_sample_id = TRUE) {

  if (length(lsets) > 1) {
    result <- purrr::map(lsets, colnames) %>% purrr::reduce(., intersect)
  } else {
    result <- intersect(colnames(lsets[[1]]), names(config))
    result <- c(result, "chemical", "endpoint")
  }
  if (remove_sample_id) {
    result <- result[!result %in% c('sample_id')]
  }


  return(result)
}

#' Nest the lsets into columns.
#'
#' @param lsets The result list from the [combi_run_rcurvep()] or [run_rcurvep()].
#' @param base_cols A vector of common column names in sets from [get_base_cols()].
#'
#' @return A tibble with the common column from sets + the data of sets nested into columns.
#' @keywords internal
#' @noRd

get_nested_joined_sets <- function(lsets, base_cols) {


  # nest the data
  lsets_n <- purrr::map2(
    lsets, names(lsets),
    function(x, y, nest_cols)
      tidyr::nest(x, !!y := -c(!!!nest_cols)), nest_cols = base_cols)

  # join all the sets
  result <- purrr::reduce(lsets_n, dplyr::left_join, by = base_cols)
  return(result)
}


#' Unnest the nested tibble with columns act_set, resp_set, fp_set, or act_summary.
#'
#'
#' @param nested A nested tibble with columns.
#' @inheritParams get_nested_joined_sets
#' @param add_col One of the column names: act_set, resp_set, fp_set, or act_summary.
#'
#' @return A tibble with unnested data.
#' @keywords internal
#' @noRd
#'
unnest_joined_sets <- function(nested, base_cols, add_col) {

  result <- nested %>%
    dplyr::select(base_cols, add_col) %>%
    tidyr::unnest(cols = add_col)
  return(result)
}


#' A wrapper function for the `add_hit_actset()`.
#'
#' @param nestd The nested tibble from the [get_nested_joined_sets()].
#'
#' @return The nested tibble but in each table of the act_set, a new column, hit is added.
#' @keywords internal
#' @noRd

add_hit <- function(nestd) {

  result <- nestd %>%
    dplyr::mutate(
      act_set = purrr::map(.data$act_set, add_hit_actset)
    )
  return(result)
}

#' Add hit column for the act_set.
#'
#' wAUC != 0 as active otherwise inactive.
#'
#' @param act_set One act_set.
#'
#' @return A tibble of the act_set with a new column, hit; 1 as active.
#' @keywords internal
#' @noRd
#'
add_hit_actset <- function(act_set) {

  result <- act_set %>%
    dplyr::mutate(
      hit = 0,
      hit = replace(.data$hit, .data$wAUC != 0, 1)
    )
  return(result)
}

#' Use the Comment and hit column to adjust the activity call.
#'
#' @param nestd The nested tibble from the [get_nested_joined_sets()]
#'   and the act_set has the hit column.
#' @inheritParams summarize_rcurvep_output
#'
#' @return The tibble but column values could be modified if inactivate is not NULL.
#' @keywords internal
#' @noRd
#'
apply_comment_to_unhit <- function(nestd, inactivate) {

  result <- nestd
  if (!is.null(inactivate)) {
    # result <- nestd %>%
    #   dplyr::mutate(
    #     act_set = purrr::map(.data$act_set, apply_comment_to_unhit_actset, inactivate = inactivate)
    #   )
    if (rlang::has_name(result, "resp_set")) {
      result <- result %>% dplyr::mutate(
        resp_set = purrr::map2(.data$resp_set, .data$act_set, apply_comment_to_unhit_nonactset))
    }
    if (rlang::has_name(result, "fp_set")) {
      result <- result %>% dplyr::mutate(
        fp_set = purrr::map2(.data$fp_set, .data$act_set, apply_comment_to_unhit_nonactset))
    }
  }
  return(result)
}


#' Use the comment to unhit values.
#'
#' @param act_set One act_set.
#' @inheritParams summarize_rcurvep_output
#'
#' @return A tibble with all activity related columns modified.
#' @keywords internal
#' @noRd
#'
apply_comment_to_unhit_actset <- function(act_set, inactivate) {
  v_cols <- c('Emax', 'slope', 'wAUC', 'wAUC_prev', 'wResp', 'AUC', 'hit')
  c_cols <- c('wConc', 'EC50', 'C50', 'POD')

  result <- act_set %>%
    dplyr::mutate_at(
      dplyr::vars(tidyselect::one_of(v_cols)),
      replace_active_value,
      inactivate = inactivate, comment = act_set$Comments, new_value = 0
    ) %>%
    dplyr::mutate_at(
      dplyr::vars(tidyselect::one_of(c_cols)),
      replace_active_value,
      inactivate = inactivate, comment = act_set$Comments, new_value = as.numeric(NA)
    )

  return(result)
}

#' Add a noted flag on the Comments column
#'
#' @param act_set a tibble of the act_set
#' @inheritParams summarize_rcurvep_output
#'
#' @return act_set with modified Commments column
#' @keywords internal
#' @noRd

adjust_rcurvep_comment <- function(act_set, inactivate) {

  if (is.character(inactivate)) {
    inp <- stringr::str_detect(act_set$Comments, inactivate)
  } else {
    inp <- as.logical(replace(rep(0, nrow(act_set)), inactivate, 1))
  }
  result <- act_set %>%
    dplyr::mutate(
      Comments = dplyr::case_when(
        inp ~ stringr::str_c(.data$Comments, "|custom"),
        TRUE ~ Comments
      )
    )
  return(result)

}

#' Replace the original value if there is a string in the comment.
#'
#' @inheritParams summarize_rcurvep_output
#' @param comment The comment string from Curvep.
#' @param ori_value The original values.
#' @param new_value The new value.
#'
#' @return A vector with some original values being modified.
#' @keywords internal
#' @noRd
#'

replace_active_value <- function(inactivate, comment, ori_value, new_value) {
  inp <- inactivate
  if (is.character(inactivate)) {
    inp <- stringr::str_detect(comment, inactivate)
  }

  result <- ori_value
  result <- replace(result, inp, new_value)
  return(result)
}


#' Use the hit column from act_set to adjust the corrected_resp in resp_set or Cxx/ECxx in fp_set.
#'
#' @param resp_set Either one resp_set or one fp_set.
#' @param act_set One act_set.
#'
#' @return The modified nonactset.
#' @keywords internal
#' @noRd
#'
apply_comment_to_unhit_nonactset <- function(nonact_set, act_set) {
  result <- nonact_set
  type <- "fp_set"
  if (rlang::has_name(result, "corrected_resp")) {
    type <- "resp_set"
  }

  if (rlang::has_name(result, "sample_id")) {

    # get the sample_ids that hit = 0
    ids <- act_set %>% dplyr::filter(.data$hit == 0) %>% dplyr::pull(.data$sample_id)

    # modify the tibble
    if (length(ids) != 0) {
      if (type == "resp_set") { #resp_set
        result <- result %>%
          dplyr::mutate(
            corrected_resp = dplyr::case_when(
              sample_id %in% ids ~ 0,
              TRUE ~ corrected_resp
            )
          )
      } else if (type == "fp_set") { #fp_set
        result <- result %>%
          dplyr::mutate(
            ECxx = dplyr::case_when(
              sample_id %in% ids ~ as.numeric(NA),
              TRUE ~ ECxx
            ),
            Cxx = dplyr::case_when(
              sample_id %in% ids ~ as.numeric(NA),
              TRUE ~ Cxx
            )
          )
      }

    }
  } else { # this is the case when there is no sample_id, only one line
    one_hit <- act_set$hit
    if (one_hit == 0) {
      if (type == "resp_set") { #resp_set
        result <- result %>% dplyr::mutate(corrected_resp = 0)
      } else if (type == "fp_set") { #fp_set
        result <- result %>% dplyr::mutate(
          ECxx = as.numeric(NA),
          Cxx = as.numeric(NA)
        )
      }

    }
  }
  return(result)
}


#' Replace all NA in the activity related columns as highest tested concentration.
#'
#' @param act_set One act_set.
#'
#' @return An act_set with modified columns (wConc, EC50, C50, POD, ECxx).
#' @keywords internal
#' @noRd

make_act_na_highconc_in <- function(act_set) {

  potency_cols <- c("wConc", "EC50", "C50", "POD", "ECxx")
  highest_conc <- act_set$highest_conc
  suppressWarnings(result <- act_set %>%
    dplyr::mutate_at(
      dplyr::vars(tidyselect::one_of(potency_cols)),
      ~ replace(.x, is.na(.x), highest_conc[is.na(.x)])
    ))

  return(result)
}




#' Replace all NA in the activity related columns as highest tested concentration (wrapper).
#'
#' @param nestd The nested tibble from the [get_nested_joined_sets()].
#'
#' @return The same nested tibble but act_set has been modified.
#' @keywords internal
#' @noRd

make_act_na_highconc <- function(nestd) {
  result <- nestd %>%
    dplyr::mutate(
      act_set = purrr::map(.data$act_set, make_act_na_highconc_in)
    )
  return(result)
}

#' Summarize act_sets
#'
#' @param nestd The nested tibble from the [get_nested_joined_sets()]
#' @param ci_level The confidence level.
#'
#' @return The original tibble with a new column, act_summary.
#' @keywords internal
#' @noRd

summarize_actsets <- function(nestd, ci_level) {
  result <- nestd %>%
    dplyr::mutate(
      act_summary = purrr::map(.data$act_set, summarize_actset_in, ci_level = ci_level)
    )

  return(result)
}

#' Summarize an active_set.
#'
#' @param act_set one active_set
#' @param ci_level confidence interval
#'
#' @return A tibble with summarized columns.
#' @keywords internal
#' @noRd
#'
summarize_actset_in <- function(act_set, ci_level) {

  # group data
  g_cols <- c('lowest_conc', 'highest_conc', 'n_conc', 'mean_conc_spacing')
  act_g_cols <- intersect(colnames(act_set), g_cols)
  act_setg <- act_set %>%
    dplyr::group_by_at(
      dplyr::vars(
        tidyselect::one_of(act_g_cols)
      )
    )

  # summarize by different types
  result <- purrr::reduce(
    list(
      summarize_actset_by_type(act_setg, type = "med_only", ci_level = ci_level),
      summarize_actset_by_type(act_setg, type = "med_conf", ci_level = ci_level),
      summarize_actset_by_type(act_setg, type = "hit_conf", ci_level = ci_level)
    ), .f = dplyr::left_join, by = act_g_cols)

  return(result)
}

#' Summarize the activity columns based on type
#'
#' @param act_set_grouped a grouped act_set
#' @param type either med_only, med_conf, hit_conf
#' @param ci_level confidence interval level
#'
#' @return a tibble with added columns (_med, _ciu, _cil, hit_confidence) depending on the type
#' @keywords internal
#' @noRd
#'
summarize_actset_by_type <- function(act_set_grouped, type, ci_level) {


  if (type == "med_only") {
    med_cols <- c('max_resp', 'min_resp', 'nCorrected')
    suppressWarnings(result <- act_set_grouped %>%
      dplyr::summarise_at(
        dplyr::vars(tidyselect::one_of(med_cols)),
        list(med = median)
      ))
  } else if (type == "med_conf") {
    upper_bound <- 1 - (1 - ci_level)/2
    lower_bound <- (1 - ci_level)/2
    ci_cols <- c('Emax', 'slope', 'AUC', 'wAUC', 'wAUC_prev','EC50','POD', 'ECxx')
    suppressWarnings(result <- act_set_grouped %>%
      dplyr::summarise_at(
        dplyr::vars(tidyselect::one_of(ci_cols)),
        list(
          med = median,
          ciu = ~ quantile(., probs = upper_bound),
          cil = ~ quantile(., probs = lower_bound)
        )
      ))
  } else if (type == "hit_conf") {
    result <- act_set_grouped %>%
      dplyr::summarize(
        n_curves = dplyr::n(),
        hit_confidence = sum(.data$hit)/dplyr::n()
      )
  }

  result <- result %>% dplyr::ungroup()
  return(result)
}
