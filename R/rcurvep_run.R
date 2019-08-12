#' curvep default
#' @export
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
    TrustHi = FALSE,
    StrictImp = TRUE
  )

  class(defaults) <- "curvep_config"
  return(defaults)
}


#'  directionality range
#'
#' @param d
#'
#' @return
#' @keywords internal
#'
#' @examples
set_curvep_range <- function(direct) {

  range <- ifelse(direct == 1, 1000000, -1000000)
  return(range)
}

create_mask <- function(d, mask = NULL) {
  if(is.null(mask)) {
    result <- d %>%
      dplyr::mutate(
        mask = 0
      )
  } else {
    result <- d %>%
      dplyr::arrange(endpoint, chemical, desc(conc)) %>%
      tidyr::nest(-endpoint, -chemical, .key = "data") %>%
      dplyr::mutate(
        mask = purrr::map(data, ~ replace(rep(0, nrow(.x)), mask, 1))
      ) %>%
      tidyr::unest()
  }
  return(result)
}

#' run  rcurvep
#'
#' @param d
#' @param config
#' @param ...
#'
#' @return
#' @export
#'
#' @examples
run_rcurvep <- function(d, mask = NULL, config = curvep_defaults(), ...) {

  d <- na.omit(d)
  d <- .check_dat_base(d)
  config <- .check_config_name(config = config, ...)
  config <- .check_config_value(config)

  d <- create_mask(d, mask)

  result <- d %>%
    dplyr::arrange(endpoint, chemical, conc) %>%
    tidyr::nest(-endpoint, -chemical, .key = "data") %>%
    dplyr::mutate(
      output = purrr::map(data, ~ call_curvep(.x$conc, .x$resp, .x$mask, config)),
      activity = purrr::map(output, tabulate_curvep_output)
    )

  return(list(result = result, config = config))
}

