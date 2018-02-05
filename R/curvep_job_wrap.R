select_type_simulation <- function(dats, directionality, n_sample, vehicle_data)
{

  #decide it's percent dataset
  if (sum(colnames(dats) %in% c("duid", "concs", "n_in", "N")) == 4)
  {
    #decide the direction of simulated resp for the percentage data
    direction <- directionality
    if (directionality == 0) direction <- 1

    #percentage data simulation
    dats <- simulate_percent_dataset(dats, direction, n_sample)

    # decide it's concs resps format
  } else if (sum(colnames(dats) %in% c("duid", "concs", "resps")) == 3)
  {
    # decide whether there is vehicle data
    if (is.null(vehicle_data))
    {
      if (!is.null(n_sample))
      {
        dats <- simulate_resp_dataset(dats, n_sample)
      } else
      {
        # add median
        dats <- dats %>%
          dplyr::group_by(duid, concs) %>%
          dplyr::summarise(resps = median(resps, na.rm = TRUE)) %>%
          dplyr::ungroup() %>%
          dplyr::mutate(repeat_id = 1)
      }

    } else
    {
      dats <- simulate_lmresp_dataset(dats, n_sample, vehicle_data)
    }

  } else {

    rlang::abort("`dats` is not an allowed dataset")
  }
  return(dats)
}


#' A wrap function to run Curvep based on types of datasets
#'
#' Given a type of datasets and parameters,
#' the function generates input for `curvep()`,
#' performs calculations, generates output, and extracts activity.
#'
#' @param dats datasets
#' @seealso \code{\link{zfishdev}}
#' @seealso \code{\link{zfishbeh}}
#'
#' @param threshold a numeric value for the presumed noise threshold or a numeric vector for threshold finding
#' @param other_paras (optional) a list of other Curvep parameters to pass on
#' @param vehicle_data NULL or a numeric vector of responses in vehicle control wells
#' @return
#' \itemize{
#'   \item dduid: dataset unique ID, it is the combination of endpoint, chemical, and directionality (from input) as well as directionality used (from calculation).
#'   \item input: a tibble, including the information of concs, resps, and parameters
#'   \item output: a list, all results from `curvep()`
#'   \item activity: a tibble, extracted activity information from output
#'   \item repeat_id: repeat id, NA for the calculation using orignal responses instead of bootstrap samples.
#'   \item thres: threshold used in the calculation
#' }
#' @seealso \code{\link{curvep}} for available Curvep parameters
#' @export
#' @importFrom magrittr "%>%"
#' @examples
#'
#' # a percentage type dataset
#'
#' data(zfishdev)
#'
#'
run_curvep_job <- function(dats, directionality = c(1, 0, -1), n_sample, threshold, other_paras = list(), vehicle_data = NULL)
{
  #arguments check
  dats <- .check_dats(dats)
  directionality <- .check_directionality(directionality)
  n_sample <- .check_n_sample(n_sample)
  threshold <- .check_threshold(threshold, directionality)
  other_paras <- .check_other_paras(other_paras)


  #create the duid column
  dats <- dats %>%
    dplyr::mutate(directionality = directionality) %>%
    tidyr::unite(duid, endpoint, chemical, directionality, sep = "#") %>%
    dplyr::arrange(duid, concs)

  # select the simulation type
  dats <- select_type_simulation(dats, directionality, n_sample, vehicle_data)

  #add a new column directionality_u
  dats <- fill_directions(dats, directionality)

  #add the mask if mask does not exist
  if (!rlang::has_name(dats, "mask")) dats$mask <- NA

  #threshold could be a list
  if (rlang::is_list(threshold))
  {
    #make list as a tibble
    thres_d <- tibble::tibble(
      directionality_u = as.integer(names(threshold)), thres = unlist(threshold)
    )
    #create a new id (dduid with direction)
    dats2 <- dats %>%
      dplyr::left_join(thres_d, by = "directionality_u") %>%
      tidyr::unite(dduid, duid, directionality_u, sep = "@", remove = FALSE) %>%
      dplyr::select(-duid)

    #prefer to be used as a dataset
    dats_si <- dats2 %>%
      split(.$repeat_id) %>%
        purrr::map(function(y) {
          result <- y %>%
            split(.$dduid) %>%
            purrr::map_df(function(x)
              create_curvep_input(
                concs = x$concs, resps = x$resps, directionality = unique(x$directionality_u),
                thres = unique(x$thres), mask = x$mask, other_paras = other_paras), .id = "dduid") %>%
            tidyr::nest(-dduid, .key = "input")
          return(result)
        })
  } else {
    # run by certain threshold
    dats2 <- dats %>%
      tidyr::unite(dduid, duid, directionality_u, sep = "@", remove = FALSE) %>%
      dplyr::select(-duid)
    thr_range <- threshold %>% rlang::set_names(.)
    dats_si <- dats2 %>%
      split(.$repeat_id) %>%
      purrr::map(function(z) {
        thr_range %>%
          purrr::map_df(function(y) {
            result <- z %>%
              split(.$dduid) %>%
              purrr::map_df(function(x) create_curvep_input(
                concs = x$concs, resps = x$resps, directionality = unique(x$directionality_u),
                thres = y, mask = x$mask, other_paras = other_paras), .id = "dduid") %>%
              tidyr::nest(-dduid, .key = "input")
            return(result)
          }, .id = "thres") %>% dplyr::mutate(thres = as.numeric(thres))
      })
  }

  dats_out <- dats_si %>%
    purrr::map(function(x) x %>% dplyr::mutate(output = purrr::map(input, run_curvep)))
  dats_out <- dats_out %>%
    purrr::map_df(function(x) x %>% dplyr::mutate(activity = purrr::map(output, tabulate_curvep_output)),.id = "repeat_id") %>%
    dplyr::mutate(repeat_id = as.numeric(repeat_id))

  if (is.null(n_sample)) dats_out <- dats_out %>% dplyr::mutate(repeat_id = NA)
  return(dats_out)

}

