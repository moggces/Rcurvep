
#' Create concentration-resonse datasets that can be applied in the `run_rcurvep()`
#'
#' The input dataset is created either by summarizing the response data
#' or by simulating the response data.
#'
#' Curvep requires 1-to-1 concentration response relationship.
#' For the dataset that does not meet the requirement, the following strategies are applied:
#'
#' ## Summary (when n_samples = NULL)
#' * For dichotomous responses, percentage is reported (n_in/N*100).
#' * For continuous responses, median value of responses per concentration is reported.
#'
#' ## Simulation (when n_samples is a positive integer)
#' * For dichotomous responses, bootstrap approach is used on the "n_in" vector
#' to create a vector of percent response.
#' * For continuous responses, options are a) direct sampling;
#' b) responses from the linear fiting using the original data + error of responses based on the supplied vehicle control data
#'
#' @param d Datasets with concentration-response data.
#'   Examples are [zfishbeh] and [zfishdev].
#' @param n_samples NULL (default) for not to simulate responses or an integer number to indicate the number of responses per concentration to simulate.
#' @param vdata NULL (default) for not to simulate responses or a vector of numeric responses in vehicle control wells to use as error.
#'   This parameter only works when n_samples is not NULL; an experimental feature.
#'
#' @return The original dataset with a new column, sample_id (if n_samples is not NULL) or
#'   the summarized dataset with columns as [zfishbeh].
#' @export
#' @seealso [run_rcurvep()]
#'
#' @examples
#'
#' # datasets with continuous response data
#' data(zfishbeh)
#'
#' ## default
#' d <- create_dataset(zfishbeh)
#'
#' ## add samples
#' d <- create_dataset(zfishbeh, n_samples = 3)
#'
#' ## add samples and vdata
#' d <- create_dataset(zfishbeh, n_samples = 3, vdata = rnorm(100))
#'
#' # dataset with dichotomous response data
#' data(zfishdev)
#'
#' ## default
#' d <- create_dataset(zfishdev)
#'
#' ## add samples
#' d <- create_dataset(zfishdev, n_samples = 3)
#'
create_dataset <- function(d, n_samples = NULL, vdata = NULL) {

  # NA is not allowed
  d <- na.omit(d)

  # make sure it is either dichotomous/continuous
  d <- .check_dat(d)
  dat_type <- assign_dataset_type(d)

  n_samples <- .check_n_samples(n_samples)
  vdata <- .check_vdata(vdata, dat_type)

  result <- create_resps(d, dataset_type = dat_type, n_samples = n_samples, vdata = vdata)

  return(result)
}

#' A wrap function to determine whether to use summarize_resps or simulate_resps.
#'
#' @inheritParams create_dataset
#' @param dataset_type Either dichotomous or continuous type.
#'
#' @return The original dataset with a new column, sample_id (if n_samples is not NULL or
#'   the summarized dataset with columns as [zfishbeh].
#' @keywords internal
#' @noRd

create_resps <- function(d, dataset_type, n_samples, vdata) {
  if (is.null(n_samples)) {
    result <- summarize_resps(d, dataset_type = dataset_type)
  } else {
    result <- simulate_resps(d, dataset_type = dataset_type, n_samples = n_samples, vdata = vdata)
  }
  return(result)
}

#' Summarize the responses.
#'
#'
#' @inheritParams create_resps
#' @return The tibble with the endpoint, chemical, conc, resp columns.
#' @keywords internal
#' @noRd

summarize_resps <- function(d, dataset_type) {

  # determine the base for group_by
  baseq <- rlang::quos(.data$endpoint, .data$chemical, .data$conc)
  if (rlang::has_name(d, "mask")) {
    baseq <- rlang::quos(.data$endpoint, .data$chemical, .data$conc, .data$mask)
  }
  d <- d %>% dplyr::group_by(!!!baseq)

  # calculation
  if (dataset_type == 'continuous') {
    result <- d %>% dplyr::summarise(resp = median(.data$resp))
  } else if (dataset_type == 'dichotomous') {
    result <- d %>% dplyr::summarise(resp = (.data$n_in/.data$N)*100)
  }

  result <- result %>% dplyr::ungroup()
  return(result)
}

#' Simulate responses.
#'
#'
#' @param d A tibble with endpoint, chemical, conc, resp
#' @inheritParams create_resps
#' @return d with sample_id and new resp columns.
#' @keywords internal
#' @noRd
#'
simulate_resps <- function(d, dataset_type, n_samples, vdata) {

  # calculate the predicted resp if needed
  if (!is.null(vdata)) d <- cal_lm_pred_resp(d)

  # determine the base for nesting
  baseq <- rlang::quos(.data$endpoint, .data$chemical, .data$conc)
  if (rlang::has_name(d, "mask")) {
    baseq <- rlang::quos(.data$endpoint, .data$chemical, .data$conc, .data$mask)
  }
  d <- d %>% tidyr::nest(data = -c(!!!baseq))


  if (dataset_type == 'continuous') {

    # sample with replacement
    if (is.null(vdata)) {
      result <- d %>%
        dplyr::mutate(
          resp = purrr::map(
            data, function(x, n_samples) sample(unlist(x$resp), size = n_samples, replace = TRUE), n_samples = n_samples
          )
        )
    } else {

      # linear fit + error
      vdata <- na.omit(replace(vdata, which(define_quantile_outliers(vdata)), NA))
      errors <- create_n_errors(vdata, n = n_samples)
      result <- d %>%
        dplyr::mutate(
          resp = purrr::map(data, function(x, n_samples) rep(unlist(x$resp), n_samples), n_samples = n_samples),
          resp = purrr::map(.data$resp, function(x, errors) x + errors, errors = errors)
        )
    }
  } else if (dataset_type == "dichotomous") {

    # bootstrap
    result <- d %>%
      dplyr::mutate(
        resp = purrr::map(
          data, function(x, n_samples) cal_percent_resps(x$n_in, x$N, times = n_samples), n_samples = n_samples
        )
      )
  }

  # clean and add sample_id
  result <- clean_add_sampleid(result, n_samples = n_samples)

  return(result)

}

#' Calculate fitted response from the linear regressinon fitting (conc x resp).
#'
#' @param d Datasets like [zfishbeh]
#' @return The tibble with an updated resp column.
#' @keywords internal
#' @noRd
#'

cal_lm_pred_resp <- function(d) {

  # determine the base
  baseq <- rlang::quos(.data$endpoint, .data$chemical)
  if (rlang::has_name(d, "mask")) {
    baseq <- rlang::quos(.data$endpoint, .data$chemical, .data$mask)
  }

  # use lm fit then predict
  result <- d %>%
    tidyr::nest(lm_input = -c(!!!baseq)) %>%
    dplyr::mutate(conc = purrr::map(.data$lm_input, ~ unique(.x$conc))) %>%
    dplyr::mutate(lm_model = purrr::map(.data$lm_input, ~ lm(resp ~ conc, data = .x))) %>%
    dplyr::mutate(
      resp = purrr::pmap(., function(...) {
        l <- list(...)
        predict(l$lm_model, data.frame(conc = l$conc))
      })
    ) %>%
    dplyr::ungroup() %>%
    dplyr::select(-.data$lm_model, -.data$lm_input) %>%
    tidyr::unnest(cols = c("conc", "resp"))
  return(result)
}

#' Assign the type of datasets.
#'
#' The column names of the datasets are used to determine the type.
#'
#' @inheritParams create_dataset
#' @return A string (either dichotomous or continuous).
#' @keywords internal
#' @noRd

assign_dataset_type <- function(d) {
  dico_cols <- c("endpoint", "chemical", "conc", "n_in", "N")
  conti_cols <- c("endpoint", "chemical", "conc", "resp")

  # check columns
  if (all(rlang::has_name(d, dico_cols))) {
    type <- 'dichotomous'
  } else if (all(rlang::has_name(d, conti_cols))) {
    type <- 'continuous'
  }
  return(type)
}



#' Calculate percent of response by bootstrap.
#'
#' @param n_in The number of incidence.
#' @param N The total number of animals.
#' @param times The bootstrap times.
#'
#' @return A vector of percent of response from the bootstrap procedure.
#' @keywords internal
#' @noRd
#'
cal_percent_resps <- function(n_in, N, times) {
  vec <- rep(0, N)
  vec <- replace(vec, sample(1:length(vec), n_in, replace = FALSE), 1)
  bd <- boot::boot(vec, sample_mean, R = times)
  return(as.vector(bd$t))
}

sample_mean <- function(x, d) {
  return(mean(x[d])*100)
}

#' Define outliers by quantile and a set cutoff.
#'
#' if x fraction of outliers exist, the Tukey approach is not suitable.
#'
#' @param vec A vector of responses
#' @param degree x IQR (default = 3)
#' @param invalid_thr (default = 0.05)
#'
#' @return A logical vec with TRUE as the outliers
#' @keywords internal
#' @noRd
#'
define_quantile_outliers <- function(vec, degree = 3, invalid_thr = 0.05) {

  q25 <- quantile(vec, probs = 0.25, na.rm = TRUE)
  q75 <- quantile(vec, probs = 0.75, na.rm = TRUE)
  IQR <- q75 - q25
  up_inf <-  q75 + degree*IQR
  low_inf <- q25 - degree*IQR
  outv <- rep(FALSE, length(vec))
  outv[(vec > up_inf | vec < low_inf) & !is.na(vec) ] <- TRUE
  if (sum(outv)/length(vec) > invalid_thr) outv <-  rep(FALSE, length(vec))
  return(outv)
}

#' create_n_errors by sample with replacement.
#'
#' @param vec A vector of responses.
#' @param n Number of samples.
#'
#' @return A vector of numeric values.
#' @keywords internal
#' @noRd

create_n_errors <- function(vec, n = 100) {

  # shift mean
  mean_vec <- mean(vec)
  vec <- vec - mean_vec

  # sample
  result <- sample(vec, size = n, replace = TRUE)
  return(result)
}



#' Clean and add sample_id column.
#'
#' @param d a tibble
#' @param n_samples Number of samples.
#'
#' @return The tibble with the sample_id column.
#' @keywords internal
#' @noRd

clean_add_sampleid <- function(d, n_samples) {
  result <- d %>%
    dplyr::mutate(sample_id = list(1:n_samples)) %>%
    dplyr::select(-data) %>%
    tidyr::unnest(cols = c("resp", "sample_id")) %>%
    dplyr::arrange(.data$sample_id)
  return(result)
}


