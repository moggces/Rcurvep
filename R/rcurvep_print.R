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

print.rcurvep <- function(x, ...) {
  x <- .check_class(x, "rcurvep", "not a rcurvep object")

  n_end <- length(unique(x$result$endpoint))
  n_chem <- length(unique(x$result$chemical))
  comp_names <- stringr::str_c(names(x), collapse = ", ")
  res_col_names <- stringr::str_c(colnames(x$result), collapse = ", ")

  data_info <- stringr::str_glue(
    "rcurvep output for a dataset with {n_end} endpoint(s) and {n_chem} chemical(s)
    object components: {comp_names}
    result has columns: {res_col_names}
    "
  )

  message("\n")
  message(data_info)
  message("\n")
}
