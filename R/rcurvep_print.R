
#' @export
print.curvep_config <- function(config, ...) {
  config <- .check_class(config, "curvep_config", "not a curvep_config object")

  message("\n")
  message("curvep configuration parameters\n")

  result <- purrr::map2(names(config), config, function(x, y) {
    new_y <- y
    if (length(y) != 1) new_y <- stringr::str_c(y, collapse = " ")
    new_x <- format(stringr::str_c(x, ":"), width = 10, justify = "right")
    stringr::str_glue(
      "    {new_x}    [{new_y}]"
    )
  })
  result <- stringr::str_c(result, collapse = "\n")
  message(result)
  message("\n")
}

#' @export
print.rcurvep <- function(x, ...) {
  x <- .check_class(x, "rcurvep", "not a rcurvep object")

  n_end <- length(unique(x$result[[1]]$endpoint))
  n_chem <- length(unique(x$result[[1]]$chemical))
  comp_names <- stringr::str_c(names(x), collapse = ", ")
  result_set_names <- stringr::str_c(names(x$result), collapse = ", ")


  data_info <- stringr::str_glue(
    "{n_end} endpoint(s) and {n_chem} chemical(s)
    Components in the list: {comp_names}
    Components in the result: {result_set_names}
    "
  )

  message("\n")
  message(data_info)
  message("\n")
}
