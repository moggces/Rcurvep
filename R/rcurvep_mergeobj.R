#' Merge results from multiple rcurvep objects
#'
#' Sometimes user may want to try multiple curvep setting and pick the one that can capture the shape (wAUC != 0).
#' The highest absolute wAUC from the chemical-endpoint(-sample_id) pair will be picked.
#'
#' @param ... rcurvep objects
#'
#' @return an updated rcurvep object with config = NULL
#' @export
#'
#' @examples
#'
#' data(zfishbeh)
#'
#' # combine default + mask
#' out1 <- combi_run_rcurvep(zfishbeh, TRSH = 10)
#' out2 <- combi_run_rcurvep(zfishbeh, TRSH = 10, mask = 1)
#' m1 <- merge_rcurvep_objs(out1, out2)
#'
#' # use same set of samples to combine
#' out1 <- combi_run_rcurvep(zfishbeh, TRSH = 10, n_samples = 2, seed = 300)
#' out2 <- combi_run_rcurvep(zfishbeh, TRSH = 10, mask = 1, n_samples = 2, seed = 300)
#' m1 <- merge_rcurvep_objs(out1, out2)
#'
merge_rcurvep_objs <- function(...) {

  objs <- list(...)

  # check input
  objs <- .check_objs_type(objs)

  # similar to get_nested_joined_sets but put wAUC out
  # maybe two functions need to be combined..
  nested <- purrr::map_df(objs, ~ fold_rcurvep_result(.x$result))

  # for each pair, pick the highest absolute wAUC
  potent_nested <- pick_row_by_wAUC(nested)

  # unnest the joined sets
  names_set <- intersect(c("act_set", "resp_set", "fp_set"), names(potent_nested))
  base_cols <- setdiff(names(potent_nested), names_set)
  lsets_m <- purrr::map(names_set, unnest_joined_sets, nested = potent_nested, base_cols = base_cols ) %>%
    rlang::set_names(names_set)

  # there is extra wAUC column
  lsets_m <- purrr::imap(lsets_m, remove_wAUC_col)


  # report the results
  result <- list(
    result = lsets_m,
    config = NULL
  )

  # class
  class(result) <- c("rcurvep", class(result))

  return(result)
}

remove_wAUC_col <- function(lset, lset_name) {
  result <- lset
  if (lset_name != "act_set") {
    result <- result %>% dplyr::select(-wAUC)
  }
  return(result)
}


#' Select the row that highest absolute wAUC
#'
#' @param d a nested tibble (act_set, resp_set, fp_set)
#'
#' @return the same structure with fewer rows
#'
#' @keywords internal
#' @noRd

pick_row_by_wAUC <- function(d) {

  # grouping; warning is when sample_id is not available
  suppressWarnings(d1 <- d %>%
    dplyr::group_by_at(
      dplyr::vars(
        tidyselect::one_of(c("endpoint", "chemical", "sample_id"))
      )
    ))

  # group sort
  result <- d1 %>%
    dplyr::arrange(
      dplyr::desc(abs(wAUC)), .by_group = TRUE) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup()

  return(result)
}


#' Create nested data
#'
#' It is similar to [get_nested_joined_sets()]).
#' Maybe to combine them later.
#' For this one, wAUC is also used as nest column
#'
#' @param lsets The result list from the [combi_run_rcurvep()] or [run_rcurvep()]
#' @return a nested tibble
#'
#' @keywords internal
#' @noRd

fold_rcurvep_result <- function(lsets) {

  # need to keep the sample_id, not summary
  base_cols <- get_base_cols(lsets, remove_sample_id = FALSE)
  lset_nested <- purrr::imap(lsets, fold_rcurvep_result_in, base_cols = base_cols)
  result <- purrr::reduce(lset_nested, dplyr::left_join, by = base_cols)

  return(result)
}

#' Create nested column for one of the three types
#'
#' @param set set data
#' @param set_name act_set, resp_set, fp_set
#' @param base_cols A vector of common column names in sets from [get_base_cols()]
#'
#' @return a tibble with a new column with nested data
#' @keywords internal
#' @noRd

fold_rcurvep_result_in <- function(set, set_name, base_cols) {

  set_nameq <- rlang::ensym(set_name)

  if (set_name == "act_set") {
    result <- set %>% tidyr::nest(!!set_nameq := -c(base_cols, "wAUC"))
  } else {
    result <- set %>% tidyr::nest(!!set_nameq := -c(base_cols))
  }
  return(result)
}
