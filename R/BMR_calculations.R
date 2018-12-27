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
#' @param p1 default = NULL, or a number to manually set the first threshold index for the distance approach
#' @param p2 default = NULL, or a number to manually set last threshold index for the distance approach
#' @param plot default = TRUE, diagnostic plot
#' @param n_endpoint_page number of endpoints to be plotted per page
#'
#' @return a tibble with added columns
#' \itemize{
#'   \item pooled_variance: the pooled varience of the potency
#'   \item thresDist: the identified threshold based on distance approach
#'   \item thresDistComment: the flag to suggest whether to use the thresDist (OK, cautionary, inadvisable)
#'   \item p1, p2, distl: the index of points for the line and the calculated distance-to-line value
#'   \item curvature: the calculated curvature value
#' }
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
  threshold = "threshold", direction = "direction", potency = "POD", p1 = NULL, p2 = NULL, plot = TRUE, n_endpoint_page = 4) {

  endpoint <- rlang::sym(endpoint)
  direction <- rlang::sym(direction)
  chemical <- rlang::sym(chemical)
  threshold <- rlang::sym(threshold)
  potency <- rlang::sym(potency)

  #calculate the pooled variance for each endpoint at a threshold
  thres_vars_pre <- cal_pooled_variances_per_endpoint_direction(
    df, endpoint = endpoint, chemical = chemical, threshold = threshold, direction = direction, potency = potency
    )

  #calculate baseline noise threshold
  select_thres <- thres_vars_pre %>%
    split(., list(.[[as.character(endpoint)]], .[[as.character(direction)]]), drop = TRUE) %>%
    purrr::map_df(
      cal_exponential_inflection, xvar = as.character(threshold), yvar = "pooled_variance", p1 = p1, p2 = p2
    )

  #generate diagonistic plot
   if (plot == TRUE)
   {
     pl <- generate_diagnostic_plot(
       select_thres,
       endpoint = endpoint, direction = direction,  threshold = threshold, n_endpoint_page = n_endpoint_page
       )
     print(pl)
   }

  return(select_thres)
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


#' Calculate the inflection point on the exponential curve using approaches
#' (currently only distance-based approach is implemented)
#'
#' @param df a tibble
#' @param xvar a character column name in the `df` to be the x-axis in the exponential curve
#' @param yvar a character column name in the `df` to be the y-axis in the expoential curve
#' @param p1 default = NULL, or a number to manually set the first threshold index for the distance approach
#' @param p2 default = NULL, or a number to manually set last threshold index for the distance approach
#'
#' @return a tibble with added columns
#' \itemize{
#'   \item thresDist: the identified threshold based on distance approach
#'   \item thresDistComment: the flag to suggest whether to use the thresDist (OK, cautionary, inadvisable)
#'   \item p1, p2, distl: the index of points for the line and the calculated distance-to-line value
#'   \item curvature: the calculated curvature value
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
#' outthres2 <- outthres %>%
#' split(., list(.$endpoint, .$direction), drop = TRUE) %>%
#' map_df(cal_exponential_inflection, xvar = "threshold", yvar = "pooled_variance")
#'
cal_exponential_inflection <- function(df, xvar, yvar, p1 = NULL, p2 = NULL) {

  #columns that will be generated
  cols_gen <- c("thresDist", "thresDistComment", "dist2l", "p1", "p2")

  #remove these columns at front
  if (sum(colnames(df) %in% cols_gen) > 1) {
    warning("the old columns (thresDist, thresDistComment, dist2l, p2, p2) will be removed")
    df <- df %>%
      dplyr::select(-dplyr::one_of(cols_gen))
  }

  #calculate the distance and the curvature
  select_thres <- df %>%
    tibble::rowid_to_column(.) %>%
    dplyr::inner_join(
      cal_dist_curva(.[[as.character(xvar)]], .[[as.character(yvar)]], p1 = p1, p2 = p2) %>%
      tibble::rowid_to_column(.),
      by = "rowid"
    ) %>%
    dplyr::mutate(
      thresDist = pick_threshold(.[[as.character(xvar)]], .$dist2l, method = "dist")
    ) %>%
    dplyr::mutate(
      thresDistComment = comment_threshold(.$thresDist, .[[as.character(xvar)]], .[[as.character(yvar)]])
    ) %>%
    dplyr::select(-rowid)

  return(select_thres)
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
#' generate_diagnostic_plot(outthres)
#'
generate_diagnostic_plot <- function(
  df, endpoint = "endpoint", direction = "direction", threshold = "threshold", n_endpoint_page = 4 ) {

  endpoint <- rlang::sym(endpoint)
  direction <- rlang::sym(direction)
  threshold <- rlang::sym(threshold)

  #make the long format
  select_thres_g <- df %>%
    dplyr::select(!!endpoint, !!direction, !!threshold, "pooled_variance", "dist2l") %>%
    tidyr::gather(type, value, c("pooled_variance", "dist2l")) %>%
    dplyr::mutate(type = ordered(type, levels =  c("pooled_variance", "dist2l") )) %>%
    dplyr::mutate(endpoint = !!endpoint, direction = !!direction, threshold = !!threshold) %>%
    dplyr::mutate(direction = as.factor(direction))

  #count how many pages
  d <- select_thres_g %>% dplyr::pull(endpoint) %>% unique
  if (length(d) < n_endpoint_page) n_endpoint_page <- length(d)
  d <- split(d, ceiling(seq_along(d)/n_endpoint_page))

  #get the slope and intercept
  lm_fitd <- df %>%
    tidyr::nest(-!!endpoint, -!!direction, -p1, -p2) %>%
    dplyr::mutate(
      temp = purrr::pmap(., function(...) {
        l <- list(...)
        get_linear_coeff_by_p1p2(l$data, l$p1, l$p2)})
    ) %>%
    dplyr::select(-data) %>%
    tidyr::unnest() %>%
    dplyr::mutate(
      type = "pooled_variance",
      type = ordered(type, levels =  c("pooled_variance", "dist2l") )
    )

  #plotting
  pl <- purrr::map(d, function(x) {
    outthres_gp <- select_thres_g %>% dplyr::filter(endpoint %in% x)
    p <- ggplot2::ggplot(outthres_gp,
                         ggplot2::aes(x = threshold, y = value, color = direction))
    p <- p + ggplot2::geom_point(size = 2) +
      ggplot2::geom_line() +
      ggplot2::geom_abline(
        data = lm_fitd, ggplot2::aes(
          slope = slope, intercept = intercept), linetype = "dashed", color = "black"
      ) +
      ggplot2::facet_grid(type ~ endpoint, scales = "free")
    return(p)
  })

  return(pl)
}

comment_threshold <- function(thresdist_v, thres, vars, cor_thres = 0.95) {
  cor_out <- cal_exponential_linear_fit(thres, vars)

  out <- rep(as.character(NA), times = length(thresdist_v))

  #based on the analysis in the data-raw/bmr_exclude
  comment <- "OK"
  if (cor_out$exponential_cor >= cor_thres & cor_out$linear_cor >= cor_thres ) {
    comment <- "cautionary"
  } else if (cor_out$exponential_cor < cor_thres & cor_out$linear_cor < cor_thres) {
    comment <- "inadvisable"
  }

  out[as.logical(thresdist_v)] <- comment

  return(out)

}


cal_dist_curva <- function(thres, pvar, p1 = NULL, p2 = NULL) {

  try(
    if (length(thres) != length(pvar)) {
      stop("the length of two vectors do not match")
    }
  )

  if (is.null(p1)) { p1 <- which.max(pvar) }
  if (is.null(p2)) { p2 <- which.min(pvar) }


  try(
    if (p1 > p2 | p1 > length(thres) | p2 > length(thres)) {
      stop("p1 and p2 do not in the range of the threshold index")
    }
  )


  d <- data.frame(thres = thres, pvar = pvar)

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

  spl <- smooth.spline(unlist(d[,1]), unlist(d[,2]))
  spl_der1 <- predict(spl, deriv = 1)
  spl_der2 <- predict(spl, deriv = 2)

  result <- tibble::tibble(
    dist2l = dists,
    #der2 = spl_der2$y,
    #der1 = spl_der1$y,
    curvature = spl_der2$y/(((1 + (spl_der1$y)^2))^3/2),
    p1 = p1, p2 = p2)

  return(result)

}


pick_threshold <- function(thres, dist2l, curvature = NULL, method = c("dist", "distcurva")) {
  out <- rep(0, times = length(dist2l)) %>% rlang::set_names(thres)

  try(
    if (method == "distcurva" & is.null(curvature)) {
      stop("curvature vector is needed")
    }
  )
  maxd <- max(dist2l, na.rm = TRUE)
  maxdi <- which.max(dist2l)
  distv <- dist2l %>% rlang::set_names(thres)
  if (method == "dist")
  {
    out[maxdi] <- 1
  } else if (method == "distcurva") { # not useful, abandoned


    curvaturev <- curvature %>% rlang::set_names(thres)
    lv <- rep(FALSE, times = length(dist2l)) %>% rlang::set_names(thres)

    lv[(maxdi):(maxdi + 3)] <- TRUE #limit to four numbers
    topl <- intersect(order(curvaturev, decreasing = TRUE), which(lv))
    topli <- 1
    maxcurv <- topl[topli] #index in the whole vector
    # exceptions:
    if (maxcurv + 1  == which.min(curvaturev)) #going forward ~ like 121
    {
      topli <- topli + 1
      maxcurv <- topl[topli]
    } else if (curvaturev[maxcurv - 1] > curvaturev[maxcurv] & maxcurv - 2 != which.max(curvaturev) ) # going backward ~ like 126
    {
      maxcurv <- maxcurv - 1
    }
    out[maxcurv] <- 1
  }
  return(out)
}


cal_exponential_linear_fit <- function(thres, vars) {

  dd <- data.frame(x = thres, y = vars)

  mod_nls <- stats::nls(y ~ exp(a + b * x),
             data = dd, start = list(a = 0, b = 0))
  mod_lm <- stats::lm(y ~ x , data = dd)

  cor_nls <- stats::cor(dd$y, stats::fitted(mod_nls))
  cor_lm <- stats::cor(dd$y, stats::fitted(mod_lm))

  return(list("exponential_cor" = cor_nls, "linear_cor" = cor_lm))

}

get_linear_coeff_by_p1p2 <- function(dd, p1, p2) {
  dd <- dd[c(p1, p2),]
  mod_lm <- stats::lm(pooled_variance ~ threshold, data = dd)
  result <- data.frame(as.list(mod_lm$coefficients)) %>%
    rlang::set_names(c("intercept", "slope"))
  return(result)
}

