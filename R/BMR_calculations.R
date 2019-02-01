#' Get baseline noise threshold through activity data calculated from simulated curves
#'
#' Given an activity dataset after `run_curvep_batch()` and `withdraw()`,
#' the function calculates the pooled variance of potency information (i.e., POD parameter) across chemicals
#' derives the baseline noise threshold, which is the lowest threshold that variance of POD is sufficiently reduced and stabilized.
#'
#' @param df a tbl
#' @param endpoint a chr that represents the column name of endpoint in `df`
#' @param chemical a chr that represents the column name of chemical in `df`
#' @param threshold a chr that represents the column name of threshold in `df`
#' @param direction a chr that represents the column name of direction in `df`
#' @param potency a chr that represents the column name of potency in `df`
#' @param plot default = TRUE, for diagnostic plots
#' @param n_endpoint_page number of endpoints to be plotted per page
#'
#' @return an object of class 'rcurvep_thres_stats' with two named tbls:
#'
#' \describe{
#'   \item{stats}{
#'      \itemize{
#'        \item y_exp_fit, y_lm_fit: the y value from the fitting results (exponential or linear)
#'        \item dist2l_exp, dist2l_raw: the calculated distance-to-line value
#'        based on the y from exponential fit or original y
#'        \item curvature: the calculated curvature value using smooth.spline fit results
#'        \item p1_raw, p2_raw: the index of points for the line using the original y values
#'      }
#'   }
#'   \item{outcome}{
#'     \itemize{
#'        \item thresDist_raw: the identified threshold based on distance approach using the original y values
#'        \item thresDist_exp: the identified threshold based on distance approach using the y values from the exponential fitting results
#'        \item thresComment: the flag to suggest whether to use the thresDist (OK, cautionary, check)
#'        \item p1_raw, p2_raw: the index of points for the line using the original y values
#'        \item cor_exp_fit, cor_lm_fit: the Pearson's correlation between the y values from fitting (exponential or linear) and the original y values
#'     }
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
#' outthres <- get_baseline_threshold(zfishdev_act)
#'
get_baseline_threshold <- function(
  df, endpoint = "endpoint", chemical = "chemical",
  threshold = "threshold", direction = "direction", potency = "POD", plot = TRUE, n_endpoint_page = 4) {

  if (!inherits(df, "rcurvep_out")) {
    warning("df is not a class of rcurvep_out")
  }

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

  result <- list(stats = select_thres_stats, outcome = select_thres_outcome)
  class(result) <- c("rcurvep_thres_stats", class(result))

  #generate diagonistic plot
   if (plot == TRUE)
   {
     pl <- plot.rcurvep_thres_stats(
       result,
       endpoint = endpoint, direction = direction,  threshold = threshold, n_endpoint_page = n_endpoint_page
       )
     print(pl)
   }

  print.rcurvep_thres_stats(result)

  return(result)
}

cal_pooled_variances_per_endpoint_direction <- function(df, endpoint, chemical, threshold, direction, potency) {

  endpoint <- rlang::sym(endpoint)
  direction <- rlang::sym(direction)
  chemical <- rlang::sym(chemical)
  threshold <- rlang::sym(threshold)
  potency <- rlang::sym(potency)

  if (sum(is.na(df[[as.character(potency)]])) > 0) {
    if ( "conc_highest" %in% colnames(df)  ) {
      df <- df %>% dplyr::mutate(!!potency := ifelse(is.na(!!potency), conc_highest, !!potency))
      warning("Potency of inactive compounds is set to the highest tested concentration")
    } else {
      stop("NA is not allowed in the potency column")
    }
  }

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


#' Calculate the knee point on the exponential-like curve using various approaches
#'
#' Two approaches are available: original value-based and exponential fit value-based.
#' Currently original value-based approach appears to be more consistent.
#'
#' @param df a tbl
#' @param xvar A chr column name in the `df` to be the x-axis in the exponential-like curve
#' @param yvar A chr column name in the `df` to be the y-axis in the exponential-like curve
#' @param p1_raw default = NULL, or an int to manually set the first threshold index for the original value-based approach
#' @param p2_raw default = NULL, or an int to manually set the last threshold index for the original value-based approach
#'
#' @return same as \code{\link{get_baseline_threshold}}
#'
#'
#' @importFrom magrittr "%>%"
#' @importFrom rlang .data
#'
#' @export
#'
#' @examples
#' data("zfishdev_act")
#' outthres <- get_baseline_threshold(zfishdev_act)[["stats"]]
#'
#' outthres2 <- outthres %>%
#'   nest(-endpoint, -direction) %>%
#'   mutate(
#'    knee_out = purrr::map(data, function(x) cal_knee_point(x, "threshold", "pooled_variance"))
#'   )
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
      y_exp_fit =  try(fitted(cal_exponential_fit(.[[xvar_c]], .[[yvar_c]])), silent = TRUE),
      y_lm_fit = fitted(cal_linear_fit(.[[xvar_c]], .[[yvar_c]]))
    ) %>%
    dplyr::mutate_at("y_exp_fit", as.numeric) %>%
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

#' @export
print.rcurvep_thres_stats <- function(thres_out) {
  outcome <- thres_out$outcome
  purrr::map(1:nrow(outcome), function(x) {
    dd <- outcome[x, ]
    ifelse(dd$direction == 1, "increasing", "decreasing")
    cat("Threshold is:", dd$thresDist_raw, "(", dd$thresComment, ")",
        "for the endpoint:", dd$endpoint, "at:",
        ifelse(dd$direction == 1, "increasing", "decreasing"), "direction\n")
  })
  invisible(thres_out)
}


#' Plots for inspecting the quality of the derived baseline noise thresholds
#'
#' @param x the output from `get_baseline_threshold()` or a tbl similar to the stats output from `get_baseline_threshold()`
#' @param ... \code{\link{get_baseline_threshold}} for available parameters
#'
#' @return a ggplot object
#'
#' @export
#'
#' @examples
#' data("zfishdev_act")
#' outthres <- get_baseline_threshold(zfishdev_act, plot = FALSE)
#' plot(outthres)
#'
plot.rcurvep_thres_stats <- function(x, ...) {

  df <- .check_plot_firstin(x)
  dots_out <- .check_plot_paras(...)
  endpoint <- dots_out$endpoint
  direction <- dots_out$direction
  threshold <- dots_out$threshold
  n_endpoint_page <- dots_out$n_endpoint_page

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

  #get the slope and intercept for raw data
  lm_fitd <- df %>%
    tidyr::nest(-endpoint, -p1_raw, -p2_raw) %>%
    dplyr::mutate(
      temp = purrr::pmap(., function(...) {
        l <- list(...)
        get_linear_coeff_by_p1p2(l$data, "threshold", "pooled_variance", l$p1_raw, l$p2_raw)})
    ) %>%
    dplyr::select(-data) %>%
    tidyr::unnest() %>%
    dplyr::mutate(
      type = "pooled_variance",
      type = ordered(type, levels =  c("pooled_variance", "dist2l_raw", "dist2l_exp") )
    )

  #get the slope and intercept for the exp fit data
  lm_fitd_exp <- df %>%
    tidyr::nest(-endpoint, -p1_raw, -p2_raw) %>%
    dplyr::mutate(
      temp = purrr::pmap(., function(...) {
        l <- list(...)
        get_linear_coeff_by_p1p2(l$data, "threshold", "y_exp_fit", 1, nrow(l$data))})
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
      ggplot2::geom_abline(
        data = lm_fitd_exp, ggplot2::aes(
          slope = slope, intercept = intercept), linetype = "dashed", color = "red"
      ) +
      ggplot2::geom_line(
        data = exp_fitd, ggplot2::aes(x = threshold, y = y_exp_fit), linetype = "solid", color = "red") +
      ggplot2::facet_grid(type ~ endpoint, scales = "free")
    return(p)
  })

  return(pl)
}


get_linear_coeff_by_p1p2 <- function(dd, xvar, yvar, p1, p2) {
  dd <- dd[c(p1, p2),]
  xvar <- rlang::sym(xvar)
  yvar <- rlang::sym(yvar)
  forma <- rlang::new_formula(yvar, xvar)

  result <- data.frame(intercept = as.numeric(NA), slope = as.numeric(NA))
  try({
    mod_lm <- lm(forma, data = dd)
    result <- data.frame(as.list(mod_lm$coefficients)) %>%
      rlang::set_names(c("intercept", "slope"))
  })

  return(result)
}

cal_exponential_fit <- function(thres, vars) {
  dd <- data.frame(x = thres, y = vars)

  mod_nls <-  nls(y ~ SSasymp(x, yf, y0, log_alpha), data = dd)

  return(mod_nls)

}

cal_linear_fit <- function(thres, vars) {
  dd <- data.frame(x = thres, y = vars)

  mod_lm <- lm(y ~ x , data = dd)

  return(mod_lm)

}

#https://dataplatform.cloud.ibm.com/analytics/notebooks/54d79c2a-f155-40ec-93ec-ed05b58afa39/view?access_token=6d8ec910cf2a1b3901c721fcb94638563cd646fe14400fecbb76cea6aaae2fb1
cal_dist2l <- function(thres, vars, p1 = NULL, p2 = NULL) {


  if (is.null(p1)) { p1 <- which.max(vars); if (length(p1) == 0) p1 <- 0}
  if (is.null(p2)) { p2 <- which.min(vars); if (length(p2) == 0) p2 <- 0 }


  try(
    if ( p1 >= p2 | p1 > length(thres) | p2 > length(thres)) {
      warning("p1 has to be smaller than p2, p1 and p2 do not in the range of the threshold index")
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
  ind_max <- thres[which.max(dist2l)]
  if ( length(ind_max) == 0 ) ind_max <- as.numeric(NA)
  return(ind_max)
}

cal_cor_original_fitted <- function(ori_vars, fitted_vars) {
  return(cor(ori_vars, fitted_vars))
}

comment_threshold <- function(exp_cor, linear_cor) {

  if (is.na(exp_cor) | is.na(linear_cor) )
  {
    return("check")
  }

  if (abs(linear_cor - exp_cor) <= 0.2) {
    if (exp_cor >= 0.90) {
      comment <- "cautionary"
    } else {
      comment <- "check"
    }
  } else {
    if (exp_cor >= 0.95) {
      comment <- "OK"
    } else {
      comment <- "cautionary"
    }
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



