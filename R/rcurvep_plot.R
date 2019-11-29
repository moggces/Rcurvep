#' Plot BMR diagnostic curves
#'
#' @param x the rcurvep_bmr object (see\code{\link{estimate_dataset_bmr}})
#' @param n_in_page number of endpoints in a page
#'
#' @return a ggplot object
#' @export
#'
#' @examples
#' data(zfishdev_act)
#'
#' # use the highest concentration to replace value of the inactive
#' sumd <- summarize_rcurvep_output(zfishdev_act, clean_only = TRUE)
#'
#' bmr_out <- estimate_dataset_bmr(sumd, plot = FALSE)
#' plot(bmr_out)
#'
plot.rcurvep_bmr <- function(x, n_in_page = 6) {

  d <- .check_class(x, "rcurvep_bmr", "not a rcurvep_bmr object")

  # get the stats component
  statsd <- d$stats

  # will add endpoint column if there is none
  statsd <- .check_bmr_statsd(statsd)

  # create a facet_endpoint column by uniting parameter columns + endpoint
  statsd1 <- add_facet_endpoint_col(statsd = statsd)

  # get the required columns
  plotd <- create_bmrplot_based(statsd = statsd1)

  # get a list of facet_endpoints per page
  endpoints_l <- get_facet_endpoints_in_page(plotd = plotd, n_in_page = n_in_page)

  # get the intercept and slope
  linecoeff <- purrr::map_df(c("ori", "exp"), get_dataset_linecoeff, statsd = statsd1)

  # generate diagnostic plots
  p <- purrr::map(endpoints_l, plot_diagnostic, plotd = plotd, lined = linecoeff)
  return(p)

}

#' Get line coefficients
#'
#' @param statsd the output from add_facet_endpoint_col()
#' @param fit_type either ori or exp
#'
#' @return a tibble with type, fit_type, intercept, slope
#' @keywords internal
#'
get_dataset_linecoeff <- function(statsd, fit_type = c("ori", "exp")) {

  nest_cols <- rlang::syms(c("facet_endpoint", "p1_ori", "p2_ori"))
  if (fit_type == "exp") {
    nest_cols <- rlang::syms(c("facet_endpoint", "p1_exp", "p2_exp"))
  }

  # get intercept and slope
  result <- statsd %>%
    tidyr::nest(-c(!!!nest_cols), .key = "input") %>%
    dplyr::mutate(
      temp = purrr::pmap(., ~ get_p1_p2_linecoeff(..4, "TRSH", "pvar", ..2, ..3))
    ) %>%
    dplyr::select(-.data$input) %>%
    tidyr::unnest()

  # add the type
  result <- result %>%
    dplyr::mutate(
      fit_type = fit_type,
      type = "pvar",
      type = ordered(.data$type, levels =  c("pvar", "dist2l_ori", "dist2l_exp") )
    )

  return(result)
}

#' Add the facet_endpoint column
#'
#' For the parameter columns (not including TRSH) and the endpoint column,
#' if the parameter has more than 1 unique values, this parameter column will be
#'
#' @param statsd the stats component from the list of estimate_dataset_bmr()
#'
#' @return a tibble with some columns are united
#' @keywords internal
#'
#'
add_facet_endpoint_col <- function(statsd) {
  id <- which(colnames(statsd) %in% "endpoint")

  # get the column names where n_distinct is more than 1
  facet_cols <- names(which(purrr::map_dbl(statsd[,1:id], dplyr::n_distinct) > 1))
  result <- statsd %>% tidyr::unite("facet_endpoint", c("endpoint",facet_cols))
  return(result)
}

#' Create the basic tibble for the BMR diagnostic plot
#'
#' @param statsd the output from add_facet_endpoint_col()
#'
#' @return a tibble facet_endpoint, TRSH, pvar, dist2l_ori, dist2l_exp, y_exp_fit
#' @keywords internal
#'
create_bmrplot_based <- function(statsd) {

  # select columns and make them as long format
  plotd <- statsd %>%
    dplyr::select("facet_endpoint", "TRSH", "pvar", "dist2l_ori", "dist2l_exp") %>%
    tidyr::gather(key = "type", value = "value", c("pvar", "dist2l_ori", "dist2l_exp"))

  # select another dataset with exponential fitted y
  exp_fitd <- statsd %>%
    dplyr::select("facet_endpoint", "TRSH", "y_exp_fit") %>%
    dplyr::mutate(type = "pvar")

  # join the two
  result <- plotd %>%
    dplyr::left_join(exp_fitd, by = c("facet_endpoint", "TRSH", "type")) %>%
    dplyr::mutate(type = ordered(.data$type, levels =  c("pvar", "dist2l_ori", "dist2l_exp")))

  return(result)
}

#' Split the facet_endpoints into pages
#'
#' @param plotd output from create_bmrplot_based
#' @param n_in_page number of datasets in a page
#'
#' @return a list, each component has the facet_endpoints in a page
#' @keywords internal
#'
get_facet_endpoints_in_page <- function(plotd, n_in_page) {
  new_ends <- unique(plotd[['facet_endpoint']])
  if (length(new_ends) < n_in_page) n_in_page <- length(new_ends)
  result <- split(new_ends, ceiling(seq_along(new_ends)/n_in_page))
  return(result)
}


#' Calculate the intercept and slope for the line by p1 and p2
#'
#' @param dd dataset
#' @param xvar x-axis
#' @param yvar y-axis
#' @param p1 an int for the first index of line
#' @param p2 an int for the last index of line
#'
#' @return a data frame with intercept and slope columns
#' @keywords internal
#'
#'
get_p1_p2_linecoeff <- function(dd, xvar, yvar, p1, p2) {
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


#' Plot the diagnostic plot
#'
#' @param plotd plotd
#' @param lined from get_p1_p2_linecoeff()
#' @param endpoints the facet_endpoint names in a page
#'
#' @return a ggplot object
#' @keywords internal

plot_diagnostic <- function(plotd, lined, endpoints) {
  # filter data
  plotdf <- plotd %>% dplyr::filter(.data$facet_endpoint %in% endpoints)
  linddf <- lined %>% dplyr::filter(.data$facet_endpoint %in% endpoints)

  # base line + dot
  p <- ggplot2::ggplot(plotdf, ggplot2::aes(x = TRSH, y = value, color = type))
  p <- p + ggplot2::geom_point(size = 2) +
    ggplot2::scale_color_manual(
      values = c("pvar" = "black", "dist2l_ori" = "black", "dist2l_exp" = "red")
    ) +
    ggplot2::geom_line()

  # exp fit y
  p <- p + ggplot2::geom_line(
    ggplot2::aes(x = TRSH, y = y_exp_fit), linetype = "solid", color = "red")

  # the p1-p2 line
  p <- p + ggplot2::geom_abline(
      data = linddf %>% dplyr::filter(.data$fit_type == "ori"), ggplot2::aes(
        slope = slope, intercept = intercept), linetype = "dashed", color = "black"
    ) +
    ggplot2::geom_abline(
      data = linddf %>% dplyr::filter(.data$fit_type == "exp"), ggplot2::aes(
        slope = slope, intercept = intercept), linetype = "dashed", color = "red"
    )

  p <- p + ggplot2::facet_grid(type ~ facet_endpoint, scales = "free")

  return(p)
}
