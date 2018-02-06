#' Unnest Curvep input, output, or activity lists
#'
#' The function unnest the input, output, or activity lists after `run_curvep_job()`
#'
#' @param c_out A tibble after `run_curvep_job()`.
#'
#' @param type A string to indicate which data to extract. Currently six types are implemented:
#'   act: all the activity related metrics such as potency and efficacy
#'   concs_hl: the highest and lowest tested concentration
#'   concs_in: a vector of input concentrations
#'   resps_out: a vector of output (clean) responses
#'   paras_in: manually input parameters
#'   paras: all the parameters used in the calculation
#'
#' @return Depending the specified type, a tibble with various columns is returned.
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
#' @seealso \code{\link{curvep}} for available Curvep parameters
#'
#'

extract_curvep_data <- function(c_out, type){


  if (type == "act")
  {
    result <- c_out %>% dplyr::select(-input, -output) %>% tidyr::unnest() %>%
      dplyr::mutate(hit = ifelse(wAUC != 0, 1, 0))
  } else if (type == "concs_hl")
  {
    result <- c_out %>%
      dplyr::select(-output, -activity) %>% tidyr::unnest() %>%
      dplyr::mutate(conc_highest = purrr::map(concs, max), conc_lowest = purrr::map(concs, min)) %>%
      dplyr::select(-concs, -resps, -paras) %>% tidyr::unnest()
  } else if (type == "concs_in" )
  {
    result <- c_out %>% dplyr::select(-output, -activity) %>% tidyr::unnest() %>%
      dplyr::select(-resps, -paras)
    return(result)
  } else if (type == "resps_out")
  {
    result <- c_out %>% dplyr::select(-input, -activity) %>%
      dplyr::mutate(resps = purrr::map(output, function(x) x$resp)) %>%
      dplyr::select(-output)
  } else if (type == "paras_in")
  {
    result <- c_out %>% dplyr::select(-output, -activity) %>% tidyr::unnest() %>%
      dplyr::select(-resps, -concs) %>%
      dplyr::mutate(temp = purrr::map(paras, function(x) x %>% purrr::flatten_df())) %>%
      dplyr::select(-paras) %>% tidyr::unnest()
  } else if (type == "paras")
  {
    result <- c_out %>% dplyr::select(-input, -activity) %>%
      dplyr::mutate(temp = purrr::map(output, ~.[['Settings']] %>% purrr::flatten_df())) %>%
      dplyr::select(-output) %>% tidyr::unnest()
  }
  return(result)
}
