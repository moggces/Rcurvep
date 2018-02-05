
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

generate_random_errors <- function(dmso, outliers, shift_mean = FALSE) {
  dmso[outliers] <- NA
  #dmso <- dmso[!is.na(dmso)]
  mean_dmso <- mean(dmso, na.rm = TRUE)
  if (shift_mean)
  {
    dmso <- dmso - mean_dmso
  }
  return(dmso)
}

# d1: incidence
# d2: total count
# direction: 1 (up); -1 (down)
boot_percent_resps <- function(d1, d2, direction = 1,  times = 100) {

  vec <- rep(0, d2)
  vec[sample(1:length(vec), d1)] <- 1
  bd <- boot::boot(vec, samplemean, R = times)
  if (direction == -1)
  {
    return((bd$t) - 100)
  } else
  {
    return((bd$t))
  }
}
samplemean <- function(x, d) {
  return(mean(x[d])*100)
}


simulate_percent_dataset <- function(dats, direction, n_sample)
{
  result <- dats %>%
    dplyr::select(duid, concs, n_in, N) %>%
    tidyr::unnest() %>%
    tidyr::nest(-duid, -concs) %>%
    dplyr::mutate(
      resps = purrr::map(data, function(x) boot_percent_resps(d1 = x$n_in, d2 = x$N, direction = direction, times = n_sample)),
      repeat_id = purrr::map(data, function(x) 1:n_sample)
    ) %>%
    dplyr::select(-data) %>%
    tidyr::unnest() %>%
    dplyr::arrange(repeat_id)
  return(result)
}

simulate_resp_dataset <- function(dats, n_sample) {

  result <- dats %>%
    dplyr::select(duid, concs, resps) %>%
    tidyr::unnest() %>%
    tidyr::nest(-duid, -concs) %>%
    dplyr::mutate(
      resps = purrr::map(data, function(x) sample(unlist(x), size = n_sample, replace = TRUE)),
      repeat_id = purrr::map(data, function(x) 1:n_sample)
    ) %>%
    dplyr::select(-data) %>%
    tidyr::unnest() %>%
    dplyr::arrange(repeat_id)
  return(result)
}

simulate_lmresp_dataset <- function(dats, n_sample, vehicle_data) {

  # clean the outliers
  dmso_clean <- generate_random_errors(
    dmso = vehicle_data,
    outliers = define_quantile_outliers(vehicle_data, degree = 3) %>% as.logical(),
    shift_mean = TRUE)

  # bootstrap the clean dmso data
  dmso_sample <- sample(dmso_clean[!is.na(dmso_clean)], size = n_sample, replace = TRUE) %>%
    as.list() %>% rlang::set_names(1:n_sample)

  # linear regression on the data
  dats <- dats %>%
    dplyr::select(duid, concs, resps) %>%
    tidyr::unnest() %>%
    tidyr::nest(-duid, .key = "lm_input") %>%
    dplyr::mutate(concs = purrr::map(lm_input, function(x) unique(x$concs))) %>%
    dplyr::mutate(lm_model = purrr::map(lm_input, function(x) {lm(resps ~ concs, data = x)})) %>%
    dplyr::rowwise() %>%
    dplyr::mutate(resps = list(predict(lm_model, data.frame(concs = concs)))) %>%
    dplyr::ungroup() %>%
    dplyr::select(-lm_model, -lm_input) %>%
    tidyr::unnest()

  # add the error to the regression data
  result <- dmso_sample %>%
    purrr::map_df(function(error) {
      result <- dats %>%
        dplyr::mutate(resps = resps + error)
      return(result)
    }, .id = "repeat_id") %>%
    dplyr::mutate(repeat_id = as.numeric(repeat_id))

  return(result)
}
