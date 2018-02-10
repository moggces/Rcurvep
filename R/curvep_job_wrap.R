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


#' A wrap function to run Curvep based on types of dataset
#'
#' Given a type of dataset and parameters,
#' the function generates input for `curvep()`,
#' performs calculations, and generates output.
#'
#' @param dats datasets such as \code{\link{zfishdev}} and \code{\link{zfishbeh}}
#' @param directionality an int value (1, 0, -1) to represent the presumed direction of responses for processing; 1 = up, -1 = down, 0 = both
#' @param n_sample NULL (using orignal data) or an int to indicate the number of curves to generate
#' @param threshold a numeric value for the presumed noise threshold or a numeric vector for threshold finding or a named (-1, 1) list to indicate the threshold of directions
#' @param other_paras a list of other Curvep parameters to pass on
#' @param vehicle_data NULL or a numeric vector of responses in vehicle control wells
#' @return
#' \itemize{
#'   \item dduid: dataset unique ID, it is the combination of endpoint, chemical, and directionality (from input) as well as directionality used (from calculation).
#'   \item input: a tibble, including the information of concs, resps, and parameters
#'   \item output: a list, all results from `curvep()`
#'   \item activity: a tibble, extracted activity information from output
#'   \item repeat_id: repeat id, NA for the calculation using orignal responses instead of simulated samples.
#'   \item thres: threshold used in the calculation
#' }
#' @seealso \code{\link{curvep}} for available Curvep parameters
#' @export
#' @importFrom magrittr "%>%"
#' @importFrom rlang .data
#' @examples
#'
#' # a percentage type dataset
#' data(zfishdev)
#' x <- zfishdev %>% split(.$endpoint)
#' outd <- run_curvep_job(x[[1]],
#'                       directionality = 1,
#'                       n_sample = 1,
#'                       threshold = 15,
#'                       other_paras = list(CARR = 20, TrustHi = TRUE))
#' # more examples are availabie
#' vignette("Rcurvep-intro")
#'
run_curvep_job <- function(dats, directionality = c(1, 0, -1), n_sample = NULL, threshold, other_paras = list(), vehicle_data = NULL)
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
    #tidyr::unite(duid, endpoint, chemical, directionality, sep = "#") %>%
    tidyr::unite(duid, endpoint, chemical, sep = "#-") %>%
    dplyr::arrange(duid, concs)

  # select the simulation type
  dats <- select_type_simulation(dats, directionality, n_sample, vehicle_data)

  #add a new column directionality_u
  dats <- fill_directions(dats, directionality)

  #add the mask if mask does not exist
  if (!rlang::has_name(dats, "mask")) dats$mask <- NA

  #create a new id (dduid with direction)
  dats2 <- dats %>%
    #dplyr::left_join(thres_d, by = "directionality_u") %>%
    tidyr::unite(dduid, duid, directionality_u, sep = "#-", remove = FALSE) %>%
    dplyr::select(-duid)

  #threshold could be a list
  if (rlang::is_list(threshold))
  {
      dats_si <- dats2 %>%
        split(.$repeat_id) %>%
        purrr::map(function(z) {
          z %>% split(.$directionality_u) %>%
            purrr::map2_df(., names(.), function(a, b) {
              thr_range <- threshold[[b]] %>% rlang::set_names(.)
              thr_range %>%
                purrr::map_df(function(y) {
                  result <- a %>%
                    split(.$dduid) %>%
                    purrr::map_df(function(x) create_curvep_input(
                      concs = x$concs, resps = x$resps, directionality = unique(x$directionality_u),
                      thres = y, mask = x$mask, other_paras = other_paras), .id = "dduid") %>%
                    tidyr::nest(-dduid, .key = "input")
                  return(result)
                }, .id = "thres") %>% dplyr::mutate(thres = as.numeric(thres))
            })
        })

  } else {

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

  #run background curvep
  dats_out <- dats_si %>%
    purrr::map(function(x) x %>% dplyr::mutate(output = purrr::map(input, run_curvep)))

  #get the activity from curvep out
  dats_out <- dats_out %>%
    purrr::map_df(function(x) x %>% dplyr::mutate(activity = purrr::map(output, tabulate_curvep_output)),.id = "repeat_id") %>%
    dplyr::mutate(repeat_id = as.numeric(repeat_id))

  #make sure repeat_id is numeric
  if (is.null(n_sample)) dats_out <- dats_out %>% dplyr::mutate(repeat_id = as.numeric(NA))

  #split the id to the original input
  dats_out <- dats_out %>% tidyr::separate(dduid, c("endpoint", "chemical", "direction"), sep = "#-")
  return(dats_out)

}

