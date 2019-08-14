
#' Create a dataset for rcurvep either by summarizing the resp data or by simulating the resp data
#'
#' Summary: (n_samples = NULL)
#' For dichotomous data, percentage is reported (n_in/N*100)
#' For continuous data, median resp per conc is reported
#'
#' Simulation: (n_samples is a positive integer)
#' For dichotomous data, bootstrap is used on the n_in vector
#' FOr continuous data, options
#'
#' @param d datasets such as \code{\link{zfishdev}} and \code{\link{zfishbeh}}
#' @param n_samples NULL (default) or an int to indicate the number of resp per conc to simulate
#' @param vdata NULL (default) or a dbl vector of responses in vehicle control wells
#'
#' @return d with sample_id (n_samples is not NULL) or d
#' @export
#'
#' @examples
#'
#'
create_dataset <- function(d, n_samples = NULL, vdata = NULL) {

  # NA is not allowed
  d <- na.omit(d)

  # make sure it is either dichotomous/continuous
  d <- .check_dat(d)
  dat_type <- assign_dataset_type(d)

  n_samples <- .check_n_samples(n_samples)
  vdata <- .check_vdata(vdata, dat_type)

  if (is.null(n_samples)) {
    result <- summarize_resps(d, dataset_type = dat_type)
  } else {
    result <- simulate_resps(d, dataset_type = dat_type, n_samples = n_samples, vdata = vdata)
  }
  return(result)
}

#' Assign dataset type
#'
#' @param d
#'
#' @return a character (dichotomous or continuous)
#' @keywords internal

assign_dataset_type <- function(d) {
  dico_cols <- c("endpoint", "chemical", "conc", "n_in", "N")
  conti_cols <- c("endpoint", "chemical", "conc", "resp")

  if (all(rlang::has_name(d, dico_cols))) {
    type <- 'dichotomous'
  } else if (all(rlang::has_name(d, conti_cols))) {
    type <- 'continuous'
  }
  return(type)
}

#' summarize_resps
#'
#' @param d
#' @param dataset_type
#'
#' @return d
#' @keywords internal

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


cal_percent_resps <- function(n_in, N, times) {
  vec <- rep(0, N)
  vec <- replace(vec, sample(1:length(vec), n_in, replace = FALSE), 1)
  bd <- boot::boot(vec, sample_mean, R = times)
  return(as.vector(bd$t))
}

sample_mean <- function(x, d) {
  return(mean(x[d])*100)
}

#' define_quantile_outliers
#'
#' if too many outliers (> 5%) no outliers are defined
#'
#' @param vec
#' @param degree
#' @param invalid_thr
#'
#' @return a logical vec with TRUE as the outliers
#' @keywords internal
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

#' create_n_errors
#'
#' @param vec
#' @param n
#'
#' @return a vector
#' @keywords internal

create_n_errors <- function(vec, n = 100) {

  # shift mean
  mean_vec <- mean(vec)
  vec <- vec - mean_vec

  # sample
  result <- sample(vec, size = n, replace = TRUE)
  return(result)
}

#' Calculate predicted resp from the linear regressinon fitting (conc x resp)
#'
#' @param d
#'
#' @return d with new resp column
#' @keywords internal

cal_lm_pred_resp <- function(d) {

  # determine the base
  baseq <- rlang::quos(.data$endpoint, .data$chemical)
  if (rlang::has_name(d, "mask")) {
    baseq <- rlang::quos(.data$endpoint, .data$chemical, .data$mask)
  }

  result <- d %>%
    tidyr::nest(-c(!!!baseq), .key = "lm_input") %>%
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
    tidyr::unnest()
  return(result)
}

#' clean_add_sampleid
#'
#' @param d
#' @param n_samples
#'
#' @return d with sample_id
#' @keywords internal

clean_add_sampleid <- function(d, n_samples) {
  result <- d %>%
    dplyr::mutate(sample_id = list(1:n_samples)) %>%
    dplyr::select(-data) %>%
    tidyr::unnest() %>%
    dplyr::arrange(.data$sample_id)
  return(result)
}

#' simulate_resps
#'
#' @param d
#' @param dataset_type
#' @param n_samples
#' @param vdata
#'
#' @return d with sample_id and new resp columns
#' @keywords internal
#'
simulate_resps <- function(d, dataset_type, n_samples, vdata = NULL) {

  # calculate the predicted resp if needed
  if (!is.null(vdata)) d <- cal_lm_pred_resp(d)

  # determine the base for nesting
  baseq <- rlang::quos(.data$endpoint, .data$chemical, .data$conc)
  if (rlang::has_name(d, "mask")) {
    baseq <- rlang::quos(.data$endpoint, .data$chemical, .data$conc, .data$mask)
  }
  d <- d %>% tidyr::nest(-c(!!!baseq), .key = "data")


  if (dataset_type == 'continuous') {

    if (is.null(vdata)) {
      result <- d %>%
        dplyr::mutate(
          resp = purrr::map(data, ~ sample(unlist(.x$resp), size = n_samples, replace = TRUE))
        )
    } else {

      vdata <- na.omit(replace(vdata, which(define_quantile_outliers(vdata)), NA))
      errors <- create_n_errors(vdata, n = n_samples)
      result <- d %>%
        dplyr::mutate(
          resp = purrr::map(data, ~ rep(unlist(.x$resp), n_samples)),
          resp = purrr::map(.data$resp, ~ .x + errors)
        )
    }
  } else if (dataset_type == "dichotomous") {

    result <- d %>%
      dplyr::mutate(
        resp = purrr::map(data, ~ cal_percent_resps(.x$n_in, .x$N, times = n_samples))
      )
  }

  # clean and add sample_id
  result <- clean_add_sampleid(result, n_samples = n_samples)

  return(result)

}
