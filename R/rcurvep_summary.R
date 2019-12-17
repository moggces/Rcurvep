#' Clean and summarize the output of rcurvep
#'
#' 1. adds an hit column in the act_set
#' 2. unhit (make result as inactive) if the Comments column contains a certain string
#' 3. make NA in the potency related columns as highest tested concentration
#' 4. summarize the results based on the ci_level
#'
#' @param d the list from the \code{\link{combi_run_rcurvep}} or \code{\link{run_rcurvep}}
#' @param act_modifier a character string, default = NULL, to make the curve with this string in the Comments column as inactive. The most common case to unhit is the "INVERSE" curves.
#' @param ci_level default = 0.95 (95 percent of confidence interval)
#' @param clean_only default = FALSE, only the 1st, 2nd, 3rd task will be performed
#'
#' @return the original list (with modified components) + act_summary
#' @export
#'
#' @examples
#'
#' data(zfishbeh)
#'
#' # base dataset
#' out <- combi_run_rcurvep(zfishbeh, n_samples = NULL, TRSH = c(5, 10))
#' out_res <- summarize_rcurvep_output(out)
#'
#' # base dataset + remove comment has "INVERSE"
#' out <- summarize_rcurvep_output(out, act_modifier = "INVERSE")
#'
#' # simulated dataset
#' out <- combi_run_rcurvep(zfishbeh, n_samples = 3, TRSH = c(5, 10))
#' out_res <- summarize_rcurvep_output(out)
#'
#'
summarize_rcurvep_output <- function(d, act_modifier = NULL, ci_level = 0.95, clean_only = FALSE) {

  d <- .check_combirun_out(d)
  lsets <- d$result
  lsets <- .check_result_sets(lsets)

  # get the common column names (e.g., TRSH)
  base_cols <- get_base_cols(lsets)

  # nest the sets into columns
  lsets_n <- get_nested_joined_sets(lsets, base_cols)

  # add the hit column in the act_set
  lsets_n <- add_hit(lsets_n)

  # use the Comments column to unhit
  lsets_n <- apply_comment_to_unhit(lsets_n, act_modifier = act_modifier)

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

#' Unnest the nested tibble with columns act_set, resp_set, fp_set, or act_summary
#'
#' @param nested a nested tibble with columns
#' @param base_cols a vector of common column names in sets
#' @param add_col act_set, resp_set, fp_set, or act_summary
#'
#' @return a tibble
#' @keywords internal
#'
unnest_joined_sets <- function(nested, base_cols, add_col) {

  result <- nested %>%
    dplyr::select(base_cols, add_col) %>%
    tidyr::unnest()
  return(result)
}

#' Get common column names between sets
#'
#' Depending on different settings (e.g., TRSH), the common column names change
#'
#' @param lsets The result list from the combi_run_rcurvep or run_rcurvep
#'
#' @return a vector of common column names in sets
#' @keywords internal

get_base_cols <- function(lsets) {

  if (length(lsets) > 1) {
    result <- purrr::map(lsets, colnames) %>% purrr::reduce(., intersect)
  } else {
    ind <- which(colnames(lsets[[1]]) %in% "chemical")
    result <- colnames(lsets[[1]])[1:ind]
  }
  result <- result[!result %in% c('sample_id')]

  return(result)
}

#' Nest the lsets into columns
#'
#' @param lsets The result list from the combi_run_rcurvep or run_rcurvep
#' @param base_cols  a vector of common column names in sets
#'
#' @return a tibble with the common column from sets + the data of sets nested into columns
#' @keywords internal

get_nested_joined_sets <- function(lsets, base_cols) {


  # nest the data
  lsets_n <- purrr::map2(
    lsets, names(lsets),
    function(x, y, nest_cols) tidyr::nest(x, -c(!!!nest_cols), .key = !!y), nest_cols = base_cols)

  # join all the sets
  result <- purrr::reduce(lsets_n, dplyr::left_join, by = base_cols)
  return(result)
}

#' A wrapper function for the add_hit_actset
#'
#' @param nestd the nested tibble from the get_nested_joined_sets
#'
#' @return the nested tibble but in each table of the act_set, a new column hit is added
#' @keywords internal

add_hit <- function(nestd) {

  result <- nestd %>%
    dplyr::mutate(
      act_set = purrr::map(.data$act_set, add_hit_actset)
    )
  return(result)
}

#' Add hit column for the act_set
#'
#' wAUC != 0 as active otherwise inactive
#'
#' @param act_set one act_set
#'
#' @return a tibble of the act_set with a new column hit; 1 as active
#' @keywords internal
#'
add_hit_actset <- function(act_set) {

  result <- act_set %>%
    dplyr::mutate(
      hit = 0,
      hit = replace(.data$hit, .data$wAUC != 0, 1)
    )
  return(result)
}

#' Use the Comment and hit column to adjust the activity call
#'
#' @param nestd the nested tibble from the get_nested_joined_sets and the act_set has the hit column
#' @param act_modifier default = NULL, no modification on the activity
#'
#' @return the same tibble but column values could be modified
#' @keywords internal
#'
apply_comment_to_unhit <- function(nestd, act_modifier) {

  result <- nestd
  if (!is.null(act_modifier)) {
    result <- nestd %>%
      dplyr::mutate(
        act_set = purrr::map(.data$act_set, apply_comment_to_unhit_actset, act_modifier = act_modifier),
        resp_set = purrr::map2(.data$resp_set, .data$act_set, apply_comment_to_unhit_nonactset),
        fp_set = purrr::map2(.data$fp_set, .data$act_set, apply_comment_to_unhit_nonactset)
      )
  }
  return(result)
}


#' Use the comment to unhit values
#'
#' @param act_set one act_set
#' @param act_modifier a character string
#'
#' @return a tibble with all activity related columns modified
#' @keywords internal
#'
apply_comment_to_unhit_actset <- function(act_set, act_modifier) {
  v_cols <- c('Emax', 'slope', 'wAUC', 'wAUC_prev', 'wResp', 'AUC', 'hit')
  c_cols <- c('wConc', 'EC50', 'C50', 'POD')

  result <- act_set %>%
    dplyr::mutate_at(
      dplyr::vars(tidyselect::one_of(v_cols)),
      replace_active_value,
      act_modifier = act_modifier, comment = act_set$Comments, new_value = 0
    ) %>%
    dplyr::mutate_at(
      dplyr::vars(tidyselect::one_of(c_cols)),
      replace_active_value,
      act_modifier = act_modifier, comment = act_set$Comments, new_value = as.numeric(NA)
    )

  return(result)
}

#' Replace the original value if there is a string in the comment
#'
#' @param act_modifier a character string
#' @param comment the comment string
#' @param ori_value original values
#' @param new_value new value
#'
#' @return  a vector with some original values being modified
#' @keywords internal

replace_active_value <- function(act_modifier, comment, ori_value, new_value) {
  result <- ori_value
  result <- replace(result, stringr::str_detect(comment, act_modifier), new_value)
  return(result)
}


#' Use the hit column from act_set to adjust the corrected_resp in resp_set or Cxx/ECxx in fp_set
#'
#' @param resp_set either one resp_set or one fp_set
#' @param act_set one act_set
#'
#' @return the modified nonactset
#' @keywords internal
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


#' Replace all NA in the activity related columns as highest tested concentration
#'
#' @param act_set one act_set
#'
#' @return an act_set with modified columns (wConc, EC50, C50, POD)
#' @keywords internal

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




#' Replace all NA in the activity related columns as highest tested concentration (wrapper)
#'
#' @param nestd the nested tibble from the get_nested_joined_sets
#'
#' @return the same nested tibble but act_set has been modified
#' @keywords internal

make_act_na_highconc <- function(nestd) {
  result <- nestd %>%
    dplyr::mutate(
      act_set = purrr::map(.data$act_set, make_act_na_highconc_in)
    )
  return(result)
}

#' Summarize act_sets
#'
#' @param nestd the nested tibble from the get_nested_joined_sets
#' @param ci_level confidence level
#'
#' @return the original tibble with a new column act_summary
#' @keywords internal

summarize_actsets <- function(nestd, ci_level) {
  result <- nestd %>%
    dplyr::mutate(
      act_summary = purrr::map(.data$act_set, summarize_actset_in, ci_level = ci_level)
    )

  return(result)
}

#' Summarize an active_set
#'
#' @param act_set one active_set
#' @param ci_level confidence interval
#'
#' @return a tibble with summarized columns
#' @keywords internal
#'
summarize_actset_in <- function(act_set, ci_level) {

  # group data
  g_cols <- c('lowest_conc', 'highest_conc', 'n_conc', 'mean_conc_spacing')
  g_cols_s <- rlang::syms(g_cols)
  act_setg <- act_set %>% dplyr::group_by(!!!g_cols_s)

  # summarize by different types
  result <- purrr::reduce(
    list(
      summarize_actset_by_type(act_setg, type = "med_only", ci_level = ci_level),
      summarize_actset_by_type(act_setg, type = "med_conf", ci_level = ci_level),
      summarize_actset_by_type(act_setg, type = "hit_conf", ci_level = ci_level)
    ), .f = dplyr::left_join, by = g_cols)

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
    ci_cols <- c('Emax', 'slope', 'AUC', 'wAUC', 'wAUC_prev','EC50','POD')
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
