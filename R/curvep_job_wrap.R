select_type_simulation <- function(dats, directionality, n_sample, vehicle_data)
{

  #decide it's percent dataset
  if (sum(colnames(dats) %in% c("duid", "concs", "n_in", "N")) == 4)
  {
    #decide the direction of simulated resp for the percentage data
    direction <- directionality
    if (directionality == 0) direction <- 1

    #percentage data simulation
    if (!is.null(n_sample))
    {
      dats <- simulate_percent_dataset(dats, direction, n_sample)
    } else
    {
      # calculate percent
      dats <- dats %>%
        dplyr::group_by(duid, concs) %>%
        dplyr::summarise(resps = (n_in/N)*100) %>%
        dplyr::mutate(
          resps = dplyr::case_when(
            direction == -1 ~ resps - 100,
            TRUE ~ resps
          )
        ) %>%
        dplyr::ungroup() %>%
        dplyr::mutate(repeat_id = 1)
    }

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


#' Run curvep in batch mode based on the type of dataset
#'
#' Given a type of dataset and parameters,
#' the function generates input for `curvep()`,
#' performs calculations, and generates output.
#'
#' \code{\link{zfishdev}} for dichotomous binary incidence data and \code{\link{zfishbeh}} for continuous data
#'
#' @param dats datasets such as \code{\link{zfishdev}} and \code{\link{zfishbeh}}
#' @param directionality an int (1, 0, -1) to represent the presumed direction of responses for processing; 1 = up, -1 = down, 0 = both
#' @param n_sample NULL (using original data) or an int to indicate the number of curves to simulate
#' @param threshold a dbl for the presumed noise threshold or a dbl vector for threshold finding or a named (-1, 1) list with dbl vectors
#' @param other_paras a list of other curvep parameters to pass on
#' @param vehicle_data NULL or a dbl vector of responses in vehicle control wells
#' @param simplify_output default = FALSE, if set TRUE (withdraw(out, "act") is performed)
#' @return if simplify_output = FALSE, an object of class 'rcurvep_out_nested' (also a tbl)
#' if simplify_output = TRUE, an object of class 'rcurvep_out' (also a tbl)
#'
#' \itemize{
#'   \item input: a tbl, including the information of concs, resps, and parameters
#'   \item output: a list, all results from `curvep()`
#'   \item activity: a tbl, activity information from output
#'   \item endpoint: input endpoint information
#'   \item chemical: input chemical information
#'   \item direction: direction used in the calculation
#'   \item threshold: threshold used in the calculation
#'   \item repeat_id: repeat id, NA for the calculation using original responses instead of the simulated responses
#' }
#' @seealso \code{\link{curvep}} for available curvep parameters
#' @seealso \code{\link{withdraw.rcurvep_out_nested}} if simplify_output = TRUE
#' @export
#' @importFrom magrittr "%>%"
#' @importFrom rlang .data
#' @examples
#'
#' # a binary incidence dataset
#' data(zfishdev)
#' x <- zfishdev %>% split(.$endpoint)
#' outd <- run_curvep_batch(x[[1]],
#'                       directionality = 1,
#'                       n_sample = 1,
#'                       threshold = 15,
#'                       other_paras = list(CARR = 20, TrustHi = TRUE))
#' # more examples are availabie
#' vignette("Rcurvep-intro")
#'
run_curvep_batch <- function(dats, directionality = c(1, 0, -1), n_sample = NULL, threshold, other_paras = list(), vehicle_data = NULL, simplify_output = FALSE)
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
    tidyr::unite(duid, endpoint, chemical, sep = "#-") %>%
    dplyr::arrange(duid, concs)

  # select the simulation type or calculate the median/percent response
  dats <- select_type_simulation(dats, directionality, n_sample, vehicle_data)

  #add a new column directionality_u
  dats <- fill_directions(dats, directionality)

  #add the mask if mask does not exist
  if (!rlang::has_name(dats, "mask")) dats$mask <- NA

  #create a new id (dduid with direction)
  dats2 <- dats %>%
    tidyr::unite(dduid, duid, directionality_u, sep = "#-", remove = FALSE) %>%
    dplyr::select(-duid)

  #threshold could be a list
  if (rlang::is_list(threshold))
  {
      dats_si <- dats2 %>%
        split(.$directionality_u) %>%
          purrr::map2_df(., names(.), function(a, b) {
          thr_range <- threshold[[b]] %>% rlang::set_names(.)
          thr_range %>%
            purrr::map_df(function(x) {
              a %>%
                tidyr::nest(-dduid, -repeat_id) %>%
                dplyr::mutate(
                  input = purrr::map(data, create_curvep_input, threshold = x, paras = other_paras)
                )
            }, .id = "threshold")
          }) %>%
        dplyr::mutate(threshold = as.numeric(threshold)) %>%
        dplyr::select(-data)


  } else {

    thr_range <- threshold %>% rlang::set_names(.)
    dats_si <- thr_range %>%
      purrr::map_df(function(x) {
        dats2 %>%
          tidyr::nest(-dduid, -repeat_id) %>%
          dplyr::mutate(
            input = purrr::map(data, create_curvep_input,
                                            threshold = x, paras = other_paras)
          )
      }, .id = "threshold") %>%
      dplyr::mutate(threshold = as.numeric(threshold)) %>%
      dplyr::select(-data)
  }

  #run background curvep
  dats_out <- dats_si %>%
                dplyr::mutate(output = purrr::map(input, run_curvep))

  #get the activity from curvep out
  dats_out <- dats_out %>%
    dplyr::mutate(activity = purrr::map(output, tabulate_curvep_output)) %>%
    dplyr::mutate(repeat_id = as.numeric(repeat_id))

  #make sure repeat_id is numeric
  if (is.null(n_sample)) dats_out <- dats_out %>% dplyr::mutate(repeat_id = as.numeric(NA))

  #split the id to the original input
  dats_out <- dats_out %>% tidyr::separate(dduid, c("endpoint", "chemical", "direction"), sep = "#-")

  #simplify the output (especially for BMR finding)
  if (simplify_output) {
    class(dats_out) <- c("rcurvep_out_nested", class(dats_out))
    dats_out <- withdraw.rcurvep_out_nested(dats_out, "act")
  } else {
    class(dats_out) <- c("rcurvep_out_nested", class(dats_out))
  }

  return(dats_out)

}

