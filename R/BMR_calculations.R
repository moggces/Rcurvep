#' Identify baseline noise threshold through activity data from simulated curves
#'
#' Provided an activity dataset after `run_curvep_job()` and `extract_curvep_data()`,
#' the function calculates the pooled variance of potency information (i.e., POD parameter) across chemicals
#' and approaches are used to derive the lowest threshold that variance of POD is sufficiently reduced and even stabilized.
#'
#' @param df a tibble
#' @param endpoint a string that represents the column name of endpoint in `df`
#' @param chemical a string that represents the column name of chemical in `df`
#' @param threshold a string that represents the column name of threshold in `df`
#' @param direction a string that represents the column name of direction in `df`
#' @param potency a string that represents the column name of potency in `df`
#' @param plot default = TRUE, diagnostic plot
#' @param n_endpoint_page number of endpoints to be plotted per page
#'
##' @return a list with two tibbles
#'
#' \describe{
#'   \item{stats}{
#'   y_exp_fit, y_lm_fit: the fitted y values (exponential or linear)
#'   dist2l_exp, dist2l_raw: the calculated distance-to-line value based on fitted y values (exponential) or raw y values
#'   curvature: the calculated curvature value using smooth.spline fit results
#'   p1_raw, p2_raw: the index of points for the line using the raw y values
#'   }
#'   \item{outcome}{
#'   thresDist_raw: the identified threshold based on distance approach using the raw y values
#'   thresDist_exp: the identified threshold based on distance approach using the y values from the fitted exponential curve
#'   thresComment: the flag to suggest whether to use the thresDist (OK, cautionary, check)
#'   p1_raw, p2_raw: the index of points for the line using the raw y values
#'   cor_exp_fit, cor_lm_fit: the Pearson's correlation between the fitted y values (exponential or linear) and the raw y values
#'   }
#' }
#'
#' @importFrom magrittr "%>%"
#' @importFrom rlang .data
#'
#' @export
#'
#' @examples
#' data("zfishdev_act")
#' acts <- zfishdev_act %>% mutate(POD = ifelse(is.na(POD), conc_highest, POD))
#' outthres <- identify_basenoise_threshold(acts)
#'
identify_basenoise_threshold <- function(
  df, endpoint = "endpoint", chemical = "chemical",
  threshold = "threshold", direction = "direction", potency = "POD", plot = TRUE, n_endpoint_page = 4) {

  endpoint <- rlang::sym(endpoint)
  direction <- rlang::sym(direction)
  chemical <- rlang::sym(chemical)
  threshold <- rlang::sym(threshold)
  potency <- rlang::sym(potency)

  #calculate the pooled variance for each endpoint at a threshold
  thres_vars_pre <- cal_pooled_variances_per_endpoint_direction(
    df, endpoint = endpoint, chemical = chemical,
    threshold = threshold, direction = direction, potency = potency
    )

  #calculate baseline noise threshold
  thres_vars_out <- thres_vars_pre %>%
    tidyr::nest(-!!endpoint, -!!direction) %>%
    dplyr::mutate(
      knee_out = purrr::map(data, function(x) cal_knee_point(x, "threshold", "pooled_variance"))
    ) %>%
    dplyr::select(-data)

  select_thres_stats <- thres_vars_out %>%
    dplyr::mutate(temp = purrr::map(knee_out, ~ .x[['stats']])) %>%
    dplyr::select(-knee_out) %>%
    tidyr::unnest()

  select_thres_outcome <- thres_vars_out %>%
    dplyr::mutate(temp = purrr::map(knee_out, ~ .x[['outcome']])) %>%
    dplyr::select(-knee_out) %>%
    tidyr::unnest()

  #generate diagonistic plot
   if (plot == TRUE)
   {
     pl <- generate_diagnostic_plot(
       select_thres_stats,
       endpoint = endpoint, direction = direction,  threshold = threshold, n_endpoint_page = n_endpoint_page
       )
     print(pl)
   }

  return(list(stats = select_thres_stats, outcome = select_thres_outcome))
}

cal_pooled_variances_per_endpoint_direction <- function(df, endpoint, chemical, threshold, direction, potency) {

  endpoint <- rlang::sym(endpoint)
  direction <- rlang::sym(direction)
  chemical <- rlang::sym(chemical)
  threshold <- rlang::sym(threshold)
  potency <- rlang::sym(potency)

  if (sum(is.na(df[[as.character(potency)]])) > 0) stop("NA is not allowed in the potency column")

  #calculate the pooled variance for each endpoint at a threshold
  result <- df %>%
    dplyr::group_by(!!endpoint, !!direction, !!chemical, !!threshold) %>%
    dplyr::summarize(sd_pod = sd(!!potency), n_rep = n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(!!endpoint, !!direction, !!threshold) %>%
    dplyr::summarize(pooled_variance  = sum((.data$sd_pod^2)*(.data$n_rep - 1))/(sum(.data$n_rep) - n())) %>%
    dplyr::ungroup()

  return(result)

}


#' Calculate the knee point on the exponential curve using various approaches
#'
#'
#'
#' @param df a tibble
#' @param xvar a character column name in the `df` to be the x-axis in the exponential curve
#' @param yvar a character column name in the `df` to be the y-axis in the exponential curve
#' @param p1_raw default = NULL, or a number to manually set the first threshold index for the distance approach using the raw y values
#' @param p2_raw default = NULL, or a number to manually set last threshold index for the distance approach the raw y values
#'
#' @return a list with two tibbles
#'
#' \describe{
#'   \item{stats}{
#'   y_exp_fit, y_lm_fit: the fitted y values (exponential or linear)
#'   dist2l_exp, dist2l_raw: the calculated distance-to-line value based on fitted y values (exponential) or raw y values
#'   curvature: the calculated curvature value using smooth.spline fit results
#'   p1_raw, p2_raw: the index of points for the line using the raw y values
#'   }
#'   \item{outcome}{
#'   thresDist_raw: the identified threshold based on distance approach using the raw y values
#'   thresDist_exp: the identified threshold based on distance approach using the y values from the fitted exponential curve
#'   thresComment: the flag to suggest whether to use the thresDist (OK, cautionary, check)
#'   p1_raw, p2_raw: the index of points for the line using the raw y values
#'   cor_exp_fit, cor_lm_fit: the Pearson's correlation between the fitted y values (exponential or linear) and the raw y values
#'   }
#' }
#'
#' @importFrom magrittr "%>%"
#' @importFrom rlang .data
#'
#' @export
#'
#' @examples
#' data("zfishdev_act")
#' acts <- zfishdev_act %>% mutate(POD = ifelse(is.na(POD), conc_highest, POD))
#' outthres <- identify_basenoise_threshold(acts)[["stats"]]
#'
#' outthres2 <- outthres %>%
#' nest(-endpoint, -direction) %>%
#' mutate(
#' knee_out = purrr::map(data, function(x) cal_knee_point(x, "threshold", "pooled_variance"))
#' )
#'

cal_knee_point <- function(df, xvar, yvar, p1_raw = NULL, p2_raw = NULL) {

  xvar_c <- as.character(xvar)
  yvar_c <- as.character(yvar)

  # get fitted y values (exponential and linear)
  # calculate distance to line based on exponential fit
  # calculate curvature
  result1 <- df %>%
    dplyr::select(dplyr::one_of(c(xvar_c, yvar_c))) %>%
    dplyr::mutate(
      y_exp_fit =  fitted(cal_exponential_fit(.[[xvar_c]], .[[yvar_c]])),
      y_lm_fit = fitted(cal_linear_fit(.[[xvar_c]], .[[yvar_c]]))
    ) %>%
    dplyr::mutate(
      dist2l_exp = cal_dist2l(.[[xvar_c]], .[['y_exp_fit']], p1 = NULL, p2 = NULL)[['dist2l']],
      curvature = cal_curvature(.[[xvar_c]], .[[yvar_c]])
    ) %>%
    tibble::rowid_to_column(.)

  # calculate distance to line based on raw values
  raw_dists <- cal_dist2l(
    df[[xvar_c]], df[[yvar_c]], p1 = p1_raw, p2 = p2_raw
  ) %>%
    dplyr::rename_all(dplyr::funs(stringr::str_c(., "_raw"))) %>%
    tibble::rowid_to_column(.)

  # join two tables
  result1 <- result1 %>%
    dplyr::inner_join(raw_dists, by = "rowid") %>%
    dplyr::select(-rowid)

  # get the thresholds based on max distance
  # calculate correlation between fitted y and original y
  # use the correlation to derive a comment on the thresholds
  result2 <- tibble::tibble(
    thresDist_raw = thres_at_max_dist2l(result1$threshold, result1$dist2l_raw),
    p1_raw = unique(result1$p1_raw),
    p2_raw = unique(result1$p2_raw),
    thresDist_exp = thres_at_max_dist2l(result1$threshold, result1$dist2l_exp),
    cor_exp_fit = cal_cor_original_fitted(result1[[yvar_c]], result1$y_exp_fit),
    cor_lm_fit = cal_cor_original_fitted(result1[[yvar_c]], result1$y_lm_fit)
  ) %>%
    dplyr::mutate(
      thresComment = comment_threshold(.$cor_exp_fit, .$cor_lm_fit)
    )

  return(list(stats = result1, outcome = result2))
}





#' Generate diagnostic plots for inspecting the quality of the identified baseline noise threshold
#'
#' @param df a tibble from the identify_basenoise_threshold() or a tibble with two required columns: pooled_variance and dist2l
#' @param endpoint a string that represents the column name of endpoint in `df`
#' @param threshold a string that represents the column name of threshold in `df`
#' @param direction a string that represents the column name of direction in `df`
#' @param n_endpoint_page number of endpoints to be plotted per page
#'
#' @return a ggplot object
#'
#' @export
#'
#' @examples
#' data("zfishdev_act")
#' acts <- zfishdev_act %>% mutate(POD = ifelse(is.na(POD), conc_highest, POD))
#' outthres <- identify_basenoise_threshold(acts, plot = FALSE)
#' generate_diagnostic_plot(outthres[['stats']])
#'
generate_diagnostic_plot <- function(
  df, endpoint = "endpoint", direction = "direction", threshold = "threshold", n_endpoint_page = 4 ) {

  endpoint <- rlang::sym(endpoint)
  direction <- rlang::sym(direction)
  threshold <- rlang::sym(threshold)

  #force the name as the same
  df <- df %>%
    dplyr::mutate(endpoint = !!endpoint, direction = !!direction, threshold = !!threshold) %>%
    dplyr::mutate(direction = as.factor(direction)) %>%
    tidyr::unite(endpoint, c("endpoint", "direction"), sep = "@")

  #make the long format
  select_thres_g <- df %>%
    dplyr::select("endpoint", "threshold", "pooled_variance", "dist2l_raw", "dist2l_exp") %>%
    tidyr::gather(type, value, c("pooled_variance", "dist2l_raw", "dist2l_exp")) %>%
    dplyr::mutate(type = ordered(type, levels =  c("pooled_variance", "dist2l_raw", "dist2l_exp") ))

  #count how many pages
  d <- select_thres_g %>% dplyr::pull(endpoint) %>% unique
  if (length(d) < n_endpoint_page) n_endpoint_page <- length(d)
  d <- split(d, ceiling(seq_along(d)/n_endpoint_page))

  #get the slope and intercept
  lm_fitd <- df %>%
    tidyr::nest(-endpoint, -p1_raw, -p2_raw) %>%
    dplyr::mutate(
      temp = purrr::pmap(., function(...) {
        l <- list(...)
        get_linear_coeff_by_p1p2(l$data, l$p1_raw, l$p2_raw)})
    ) %>%
    dplyr::select(-data) %>%
    tidyr::unnest() %>%
    dplyr::mutate(
      type = "pooled_variance",
      type = ordered(type, levels =  c("pooled_variance", "dist2l_raw", "dist2l_exp") )
    )

  #add exp_fit results
  exp_fitd <- df %>%
    dplyr::select(endpoint, threshold, y_exp_fit) %>%
    dplyr::mutate(
      type = "pooled_variance",
      type = ordered(type, levels =  c("pooled_variance", "dist2l_raw", "dist2l_exp") )
    )

  #plotting
  pl <- purrr::map(d, function(x) {
    outthres_gp <- select_thres_g %>% dplyr::filter(endpoint %in% x)
    p <- ggplot2::ggplot(outthres_gp,
                         ggplot2::aes(x = threshold, y = value, color = type))
    p <- p + ggplot2::geom_point(size = 2) +
      ggplot2::scale_color_manual(
        values = c("pooled_variance" = "black", "dist2l_raw" = "black", "dist2l_exp" = "red")
      ) +
      ggplot2::geom_line() +
      ggplot2::geom_abline(
        data = lm_fitd, ggplot2::aes(
          slope = slope, intercept = intercept), linetype = "dashed", color = "black"
      ) +
      ggplot2::geom_line(
        data = exp_fitd, ggplot2::aes(x = threshold, y = y_exp_fit), linetype = "dashed", color = "red") +
      ggplot2::facet_grid(type ~ endpoint, scales = "free")
    return(p)
  })

  return(pl)
}


get_linear_coeff_by_p1p2 <- function(dd, p1, p2) {
  dd <- dd[c(p1, p2),]
  mod_lm <- lm(pooled_variance ~ threshold, data = dd)
  result <- data.frame(as.list(mod_lm$coefficients)) %>%
    rlang::set_names(c("intercept", "slope"))
  return(result)
}

cal_exponential_fit <- function(thres, vars) {
  dd <- data.frame(x = thres, y = vars)

  mod_nls <- nls(y ~ exp(a + b * x),
                        data = dd, start = list(a = 0, b = 0))

  return(mod_nls)

}

cal_linear_fit <- function(thres, vars) {
  dd <- data.frame(x = thres, y = vars)

  mod_lm <- lm(y ~ x , data = dd)

  return(mod_lm)

}

#https://dataplatform.cloud.ibm.com/analytics/notebooks/54d79c2a-f155-40ec-93ec-ed05b58afa39/view?access_token=6d8ec910cf2a1b3901c721fcb94638563cd646fe14400fecbb76cea6aaae2fb1
cal_dist2l <- function(thres, vars, p1 = NULL, p2 = NULL) {


  if (is.null(p1)) { p1 <- which.max(vars) }
  if (is.null(p2)) { p2 <- which.min(vars) }


  try(
    if (p1 > p2 | p1 > length(thres) | p2 > length(thres)) {
      warning("p1 and p2 do not in the range of the threshold index")
      dists <- rep(as.numeric(NA), length(thres))
    }
  )


  d <- data.frame(x = thres, y = vars)


  while (p1 < p2)
  {
    l1 <- d[p1, ] %>% purrr::as_vector(.)
    l2 <- d[p2, ] %>% purrr::as_vector(.)

    dists <- apply(d, 1, function(x) {
      x <- purrr::as_vector(x)
      v1 <- l1 - l2
      v2 <- x - l1
      m <- cbind(v1,v2)
      d <- abs(det(m))/sqrt(sum(v1*v1))
      return(d)
    })

    dists[!(seq(1:length(dists)) %in% p1:p2)] <- NA
    l3 <- d[which.max(dists),] %>% purrr::as_vector(.)
    if (abs(l3[2] - l2[2]) < abs(l3[2] - l1[2])) {
      break
    } else
    {
      p1 <- p1 + 1
    }
  }

  result <- tibble::tibble(
    dist2l = dists,
    p1 = p1, p2 = p2)

  return(result)
}

thres_at_max_dist2l <- function(thres, dist2l) {
  return(thres[which.max(dist2l)])
}

cal_cor_original_fitted <- function(ori_vars, fitted_vars) {
  return(cor(ori_vars, fitted_vars))
}

comment_threshold <- function(exp_cor, linear_cor) {
  comment <- "OK"
  if (exp_cor >= 0.95 & linear_cor >= 0.95 ) {
    comment <- "cautionary"
  } else if (exp_cor < 0.95 & abs(linear_cor - exp_cor) < 0.2 ) {
    comment <- "check"
  }
  return(comment)
}


cal_curvature <- function(thres, vars) {

  spl <- smooth.spline(thres, vars)
  spl_der1 <- predict(spl, deriv = 1)
  spl_der2 <- predict(spl, deriv = 2)
  curvature <-  spl_der2$y/(((1 + (spl_der1$y)^2))^3/2)
  return(curvature)
}



