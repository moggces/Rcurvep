#' Estimate benchmark response (BMR) for each dataset
#'
#' Currently two methods have been implemented to get the "keen-point" from the variance(y) - threshold(x) curve.
#' One is to use the original y values to draw a straight line between the lowest x value (p1) to highest x value (p2).
#' The knee-point is the x that has the longest distance to the line.
#' The other one is to fit the data first then use the fitted responses to do the same analysis.
#' Currently the first method is preferred.
#'
#' @details
#'
#' The estimated BMR can be used in the calculation of POD.
#' For example, if bmr = 25.
#' For Curvep, `combi_run_rcurvep(zfishbeh, TRSH = 25)`.\cr
#' For Hill fit, `summarize_fit_output(run_fit(zfishbeh, modls = "hill"), thr_resp = 25, extract_only = TRUE)`.
#'
#' @param d The rcurvep object with multiple samples and TRSHs. See [combi_run_rcurvep()] for an example.
#' @param p1 Default = NULL, or an integer value to manually set the first index of line.
#' @param p2 Default = NULL, or an integer value to manually set the last index of line.
#' @param plot Default = TRUE, plot the diagnostic plot.
#' @return A list with two components: stats and outcome.
#'
#' \itemize{
#'   \item stats: a tibble, including pooled variance (pvar),
#'   fitted responses (y_exp_fit, y_lm_fit), distance to the line (dist2l)
#'   \item outcome: a tibble, including estimated BMRs (bmr)
#' };
#' Suffix in the **stats** and **outcome** tibble: "ori" (original values), "exp"(exponential fit).
#' prefix in the **outcome** tibble, "cor" (correlation between the fitted responses and the original responses),
#' "bmr" (benchmark response), "qc" (quality control).
#'
#' @export
#' @seealso [cal_knee_point()], [combi_run_rcurvep()]
#' @examples
#'
#' # no extra cleaning
#' data(zfishdev_act)
#' bmr_out <- estimate_dataset_bmr(zfishdev_act, plot = FALSE)
#' plot(bmr_out)
#'
#' # if want to do extra cleaning...
#' actm <- summarize_rcurvep_output(zfishdev_act, clean_only = TRUE, inactivate = "CARRY_OVER")
#'
#' bmr_out <- estimate_dataset_bmr(actm, plot = FALSE)
#'
#'
estimate_dataset_bmr <- function(d, p1 = NULL, p2 = NULL, plot = TRUE) {

  # check input data
  d <- .check_bmr_input(d)


  # get the base columns (parameter columns + endpoint + chemical, no sample_id)
  base_cols <- get_base_cols(d$result)

  # get the activity data
  act_set <- make_act_na_highconc_in(d$result$act_set) # make act NA as highest conc

  # calculate the pooled variance for each dataset
  pvar_datasets <- cal_dataset_pvar(act_set, base_cols)

  # calculate knee point for each dataset
  knees <- cal_dataset_knee(pvard = pvar_datasets, base_cols = base_cols, p1 = p1, p2 = p2)

  # unnest the knee result
  resultl <- purrr::map(c("stats", "outcome"), unnest_knee_data, kneed = knees)

  # output and set class
  result <- list(stats = resultl[[1]], outcome = resultl[[2]])
  class(result) <- c("rcurvep_bmr", class(result))

  # plot the diagnostic plot
  if (plot) {
    p <- plot.rcurvep_bmr(result)
    print(p)
  }

  return(result)

}

#' Calculate pooled variance column
#'
#' @param act_set The act_set. There
#' @param base_cols parameter columns + endpoint + chemical, no sample_id
#'
#' @return a tibble with parameter columns + endpoint + pvar
#' @keywords internal
#' @noRd

cal_dataset_pvar <- function(act_set, base_cols) {

  # quo base_cols
  base_colsq <- rlang::syms(base_cols)
  base_cols_nochemq <- rlang::syms(base_cols[!base_cols %in% "chemical"])

  # calculate SD and n for each chemical
  sdd <- act_set %>%
    dplyr::group_by(!!!base_colsq) %>%
    dplyr::summarise(
      sd_pod = sd(.data$POD), n_rep = dplyr::n()
    ) %>%
    dplyr::ungroup()

  # calculate pooled variance
  result <- sdd %>%
    dplyr::group_by(!!!base_cols_nochemq) %>%
    dplyr::summarize(
      pvar  = sum((.data$sd_pod^2)*(.data$n_rep - 1))/(sum(.data$n_rep) - dplyr::n())) %>%
    dplyr::ungroup()

  return(result)
}

#' Calculate knee point for each dataset
#'
#' @param pvard output `cal_dataset_pvar()`.`
#' @param base_cols parameter columns + endpoint + chemical, no sample_id.
#' @param p1 NULL, or an int to manually set the first index of line.
#' @param p2 NULL, or an int to manually set the last index of line.
#'
#' @return the pvard + a new column knee_out
#' @keywords internal
#' @noRd
#'
#'
cal_dataset_knee <- function(pvard, base_cols, p1, p2) {
  base_cols_f <- base_cols[!base_cols %in% c("chemical", "TRSH")]
  knees <- pvard %>%
    tidyr::nest(input = -c(base_cols_f)) %>%
    dplyr::mutate(
      knee_out = purrr::map(
        .data$input,
        cal_knee_point, xaxis = "TRSH", yaxis = "pvar", p1 = p1, p2 = p2, plot = FALSE)
    )
  return(knees)
}

#' Unnest the knee_out object
#'
#' @param kneed output from cal_dataset_knee()
#' @param type either stats or outcome
#'
#' @return a tibble
#' @keywords internal
#' @noRd
#'
unnest_knee_data <- function(kneed, type = c("stats", "outcome")) {
  result <- kneed %>%
    dplyr::mutate(temp = purrr::map(.data$knee_out, ~ .x[[type]])) %>%
    dplyr::select(-.data$knee_out, -.data$input) %>%
    tidyr::unnest(cols = c("temp"))
  return(result)
}


#' Calculate the knee point on the exponential-like curve
#'
#'
#' @param d A tibble.
#' @inheritParams estimate_dataset_bmr
#' @param xaxis The column name in the `d` to be the x-axis in the exponential-like curve
#' @param yaxis The column name in the `d` to be the y-axis in the exponential-like curve
#'
#' @inherit estimate_dataset_bmr description return
#' @export
#' @seealso [estimate_dataset_bmr()]
#'
#' @examples
#'
#' inp <- data.frame(
#' x = seq(5, 95, by = 5),
#' y = c(0.0537, 0.0281, 0.0119, 0.0109, 0.0062, 0.0043, 0.0043, 0.0042,
#' 0.0041, 0.0043, 0.0044, 0.0044, 0.0046, 0.0051,
#' 0.0055, 0.0057, 0.0072, 0.0068, 0.0035)
#' )
#'
#' out <- cal_knee_point(inp,"x", "y", plot = FALSE)
#' plot(out)
#'
#'
cal_knee_point <- function(d, xaxis, yaxis, p1 = NULL, p2 = NULL, plot = TRUE) {

  xvar_c <- as.character(xaxis)
  yvar_c <- as.character(yaxis)


  # remove all the other unrelated columns
  d1 <- d %>% dplyr::select(dplyr::one_of(c(xvar_c, yvar_c)))

  # add the y
  res_stat <- add_fitted_y(d = d1, xvar = xvar_c, yvar = yvar_c)
  res_stat <- add_dist_2_line(res_stat, xvar = xvar_c, yvar = yvar_c, p1 = p1, p2 = p2)
  res_stat <- add_curvature(res_stat, xvar = xvar_c, yvar = yvar_c)

  res_out <- create_bmr_report(res_stat, xvar = xvar_c, yvar = yvar_c)

  result <- list(stats = res_stat, outcome = res_out)
  class(result) <- c("rcurvep_bmr", class(result))
  # plot the diagnostic plot
  if (plot) {
    p <- plot.rcurvep_bmr(result)
    print(p)
  }

  return(result)
}


#' Create the output for BMR from two approaches
#'
#' @param distd the output of add_dist_2_line()
#' @param xvar a chr column name
#' @param yvar a chr column name
#'
#' @return a tibble with generate columns
#' @keywords internal
#' @noRd
#'

create_bmr_report <- function(distd, xvar, yvar) {

  result <- tibble::tibble(
    bmr_ori = get_thres_at_max_dist2l(distd[[xvar]], distd[['dist2l_ori']]),
    p1_ori = unique(distd$p1_ori),
    p2_ori = unique(distd$p2_ori),
    bmr_exp = get_thres_at_max_dist2l(distd[[xvar]],distd[['dist2l_exp']]),
    p1_exp = unique(distd$p1_exp),
    p2_exp = unique(distd$p2_exp),
    cor_exp_fit = cal_cor_ori_fitted(distd[[yvar]], distd[['y_exp_fit']]),
    cor_lm_fit = cal_cor_ori_fitted(distd[[yvar]], distd[['y_lm_fit']])
  )

  result[['qc']] <- qc_data(result[['cor_exp_fit']], result[['cor_lm_fit']])
  return(result)

}


#' Quality control on the data (pvar ~ TRSH) for getting the BMR
#'
#' By using the correlation coefficient between pvar and y_exp_fit and y_lm_fit,
#' to heuristically QC the data
#'
#' @param exp_cor correlation coefficient between pvar and y_exp_fit
#' @param linear_cor correlation coefficient between pvar and y_lm_fit
#'
#' @return a character string either of three values: OK, cautionary, check
#' @keywords internal
#' @noRd

qc_data <- function(exp_cor, linear_cor) {

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


#' Add fitted values by either exponential fit or linear fit to the original table
#'
#' @param d a tibble with two columns, xvar and yvar
#' @param xvar a chr column name
#' @param yvar a chr column name
#'
#' @return the original d + y_exp_fit + y_lm_fit columns
#' @keywords internal
#' @noRd
#'
add_fitted_y <- function(d, xvar, yvar) {
  result <- d %>%
    dplyr::mutate(
      y_exp_fit =  try(fitted(cal_exponential_fit(.data[[xvar]], .data[[yvar]])), silent = TRUE),
      y_lm_fit = fitted(cal_linear_fit(.data[[xvar]], .data[[yvar]]))
    ) %>%
    dplyr::mutate_at("y_exp_fit", as.numeric)
  return(result)
}

#' Add distance to line at each x-axis point
#'
#'
#' @param fitd a tibble from the add_fitted_y
#' @param xvar a chr column name
#' @param yvar a chr column name
#' @param p1 NULL, or an int to manually set the first index of line
#' @param p2 NULL, or an int to manually set the last index of line
#'
#' @return a original fitd + 6 new columns
#' @keywords internal
#' @noRd
#'
add_dist_2_line <- function(fitd, xvar, yvar, p1, p2) {

  # distance calculation
  exp_dist <- cal_dist2l(fitd[[xvar]], fitd[['y_exp_fit']], p1 = p1, p2 = p2)
  ori_dist <- cal_dist2l(fitd[[xvar]], fitd[[yvar]], p1 = p1, p2 = p2)

  # clean data
  result <- fitd %>%
    dplyr::mutate(
      dist2l_ori = ori_dist$dist2l,
      p1_ori = ori_dist$p1,
      p2_ori = ori_dist$p2,
      dist2l_exp = exp_dist$dist2l,
      p1_exp = exp_dist$p1,
      p2_exp = exp_dist$p2

    )
  return(result)
}

#' Add curvature
#'
#' @param d a tibble with two columns, xvar and yvar
#' @param xvar a chr column name
#' @param yvar a chr column name
#'
#' @return a tibble with a new column, curva
#' @keywords internal
#' @noRd
#'
add_curvature <- function(d, xvar, yvar) {
  result <- d %>%
    dplyr::mutate(
      curva = cal_curvature(.data[[xvar]], .data[[yvar]])
    )
  return(result)
}



#' Calculate exponential fit (y ~ x)
#'
#' http://douglas-watson.github.io/post/2018-09_exponential_curve_fitting/&sa=D&source=hangouts&ust=1547125299019000&usg=AFQjCNGzoP26jIUYQnAS9Z5Eb4aO66IzPg
#'
#' @param thres x-axis
#' @param vars y-axis
#'
#' @return an nls object
#' @keywords internal
#' @noRd
#'
cal_exponential_fit <- function(thres, vars) {

  dd <- data.frame(x = thres, y = vars)
  mod_nls <-  nls(y ~ SSasymp(x, yf, y0, log_alpha), data = dd)

  return(mod_nls)

}


#' Calcualte linear fit (y ~ x)
#'
#' @param thres x-axis
#' @param vars y-axis
#'
#' @return an lm object
#' @keywords internal
#' @noRd
#'
cal_linear_fit <- function(thres, vars) {

  dd <- data.frame(x = thres, y = vars)
  mod_lm <- lm(y ~ x , data = dd)

  return(mod_lm)

}


#' Calculate curvature
#'
#' https://en.wikipedia.org/wiki/Curvature
#'
#' @param thres x-axis
#' @param vars y-axis
#'
#' @return a vector of curvature at each x point
#' @keywords internal
#' @noRd

cal_curvature <- function(thres, vars) {

  spl <- smooth.spline(thres, vars)
  spl_der1 <- predict(spl, deriv = 1)
  spl_der2 <- predict(spl, deriv = 2)
  curvature <-  spl_der2$y/(((1 + (spl_der1$y)^2))^3/2)
  return(curvature)
}


#' Calculate the distance to a line
#'
#' https://dataplatform.cloud.ibm.com/analytics/notebooks/54d79c2a-f155-40ec-93ec-ed05b58afa39/view?access_token=6d8ec910cf2a1b3901c721fcb94638563cd646fe14400fecbb76cea6aaae2fb1
#'
#' @param thres x-axis
#' @param vars y-axis
#' @param p1 NULL, or an int to manually set the first index of line
#' @param p2 NULL, or an int to manually set the last index of line
#'
#' @return a vector of curvature at each x point
#' @keywords internal
#' @noRd

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
    l1 <- d[p1, ] %>% purrr::as_vector()
    l2 <- d[p2, ] %>% purrr::as_vector()

    dists <- apply(d, 1, function(x) {
      x <- purrr::as_vector(x)
      v1 <- l1 - l2
      v2 <- x - l1
      m <- cbind(v1,v2)
      d <- abs(det(m))/sqrt(sum(v1*v1))
      return(d)
    })

    dists[!(seq(1:length(dists)) %in% p1:p2)] <- NA
    l3 <- d[which.max(dists),] %>% purrr::as_vector()
    if (abs(l3[2] - l2[2]) < abs(l3[2] - l1[2])) {
      break
    } else
    {
      p1 <- p1 + 1
    }
  }

  result <- list(dist2l = dists, p1 = p1, p2 = p2)

  return(result)
}

#' Get the threshold at which has the max distance to the line
#'
#' @param thres a vector of thresholds
#' @param dist2l a vector of distances
#'
#' @return a single value from the thres
#' @keywords internal
#' @noRd
#'
get_thres_at_max_dist2l <- function(thres, dist2l) {
  ind_max <- thres[which.max(dist2l)]
  if ( length(ind_max) == 0 ) ind_max <- as.numeric(NA)
  return(ind_max)
}

cal_cor_ori_fitted <- function(ori_vars, fitted_vars) {
  return(cor(ori_vars, fitted_vars))
}

