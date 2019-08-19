summarize_result_sets <- function(lsets, act_modifier = NULL) {
  lsets <- .check_result_sets(lsets)

  inter_cols <- get_base_cols(lsets)
  result <- get_nested_joined_sets(lsets, base_cols = inter_cols)

  result <- add_hit(result)
  result <- apply_comment_to_unhit(result, act_modifier = act_modifier)


  return(result)
}

get_base_cols <- function(lsets) {
  result <- purrr::map(lsets, colnames) %>% purrr::reduce(., intersect)
  result <- result[!result %in% c('sample_id')]
  return(result)
}

get_nested_joined_sets <- function(lsets, base_cols) {
  base_colsq <- rlang::enquos(base_cols)
  lsets_n <- purrr::map2(lsets, names(lsets),
                         function(x, y, nest_cols) tidyr::nest(x, -c(!!!nest_cols), .key = !!y), nest_cols = base_colsq)
  result <- purrr::reduce(lsets_n, dplyr::inner_join, by = base_cols)
  return(result)
}

add_hit <- function(d) {

  result <- d %>%
    dplyr::mutate(
      act_set = purrr::map(.data$act_set, add_hit_actset)
    )
  return(result)
}

add_hit_actset <- function(act_set) {

  result <- act_set %>%
    dplyr::mutate(
      hit = ifelse(.data$wAUC != 0, 1, 0)
    )
  return(result)
}

apply_comment_to_unhit <- function(d, act_modifier = NULL) {

  result <- d
  if (!is.null(act_modifier)) {
    result <- d %>%
      dplyr::mutate(
        act_set = purrr::map(.data$act_set, apply_comment_to_unhit_actset, act_modifier = act_modifier),
        resp_set = purrr::map2(.data$resp_set, .data$act_set, apply_comment_to_unhit_respset),
        fp_set = purrr::map2(.data$fp_set, .data$act_set, apply_comment_to_unhit_fpset)
      )
  }
  return(result)
}

# act_modifier is not null
apply_comment_to_unhit_actset <- function(act_set, act_modifier) {
  v_cols <- c('Emax', 'slope', 'wAUC', 'wAUC_prev', 'wResp', 'AUC', 'hit')
  c_cols <- c('wConc', 'EC50', 'C50', 'POD')


  result <- act_set %>%
    dplyr::mutate_at(
      dplyr::vars(tidyselect::one_of(v_cols)),
        function(x) ifelse(stringr::str_detect(.$Comments, act_modifier), 0, x)
    ) %>%
    dplyr::mutate_at(
      dplyr::vars(tidyselect::one_of(c_cols)),
      function(x) ifelse(stringr::str_detect(.$Comments, act_modifier), as.numeric(NA), x)
    )

  return(result)
}

# not right,
apply_comment_to_unhit_respset <- function(resp_set, act_set) {
  result <- resp_set
  ids <- which(act_set$hit == 0)
  if (rlang::has) {
    result <- result %>%
      dplyr::mutate(
        corrected_resp = 0
      )
  }
  return(result)
}

# not right,
apply_comment_to_unhit_fpset <- function(fp_set, act_set) {
  result <- fp_set
  if (act_set$hit == 0) {
    result <- result %>%
      dplyr::mutate(
        ECxx = as.numeric(NA),
        Cxx = as.numeric(NA)
      )
  }
  return(result)
}

# make_act_na_highconc <- function(d) {
#   d %>%
#     dplyr::mutate(
#       act_set = purrr::map(act_set, )
#     )
# }
#
# add_confidence_actset <- function(act_set) {
#   u_cols <-
#
#   if (rlang::has_name(act_set, "sample_id")) {
#
#   }
# }

add_median_
