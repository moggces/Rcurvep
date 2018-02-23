#' Identify baseline noise threshold through activity data from simulated curves
#'
#' Provided an activity dataset after `run_curvep_job()` and `extract_curvep_data()`,
#' the function calculates the pooled variance of potency information (i.e., POD parameter) across chemicals
#' and approaches are used to derive the lowest threshold that variance of POD is sufficiently reduced and even stabilized.
#'
#' @param df a tibble
#' @param id a vector of strings of the column names in `df`, representing the investigated endpoint
#' @param chemical a string that represents the column name of chemical in `df`
#' @param threshold a string that represents the column name of threshold in `df`
#' @param potency a string that represents the column name of potency in `df`
#'
#' @return a tibble with added columns
#' \itemize{
#'   \item pooled_variance: the pooled varience of the potency
#'   \item thresDist: the identified threshold based on distance approach
#'   \item thresCurva: the identified threshold based on curvature + distance approach
#'   \item p1, p2, distl: the index of points for the line and the calculated distance-to-line value
#'   \item curvature: the calculated curvature value
#' }
#' @importFrom magrittr "%>%"
#' @importFrom rlang .data
#'
#' @export
#'
#' @examples
#' vignette("Rcurvep-intro")
#'
identify_basenoise_threshold <- function(df, id, chemical, threshold, potency) {


  id <- rlang::syms(id)
  chemical <- rlang::sym(chemical)
  threshold <- rlang::sym(threshold)
  potency <- rlang::sym(potency)

  thres_vars_pre <- df %>%
    dplyr::group_by(!!!id, !!chemical, !!threshold) %>%
    dplyr::summarize(sd_pod = sd(!!potency), n_rep = n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(!!!id, !!threshold) %>%
    dplyr::summarize(pooled_variance  = sum((.data$sd_pod^2)*(.data$n_rep - 1))/(sum(.data$n_rep) - n())) %>%
    dplyr::ungroup()


   thres_vars <- thres_vars_pre %>%
     tidyr::nest(-c(!!!id)) %>%
     dplyr::mutate(
       outd = purrr::map(data, function(x) cal_dist_curva(x[[as.character(threshold)]], x$pooled_variance))
     ) %>%
     tidyr::unnest()

   select_thres <- thres_vars %>%
    tidyr::nest(-c(!!!id)) %>%
    dplyr::mutate(
      thresCurva = purrr::map(data, function(x) pick_threshold(x[[as.character(threshold)]], x$dist2l, x$curvature, method = "distcurva")),
      thresDist = purrr::map(data, function(x) pick_threshold(x[[as.character(threshold)]], x$dist2l, method = "dist"))

    ) %>%
    tidyr::unnest()

  return(select_thres)
}


cal_dist_curva <- function(thres, pvar) {

  try(
    if (length(thres) != length(pvar)) {
      stop("the length of two vectors do not match")
    }
  )
  p1 <- which.max(pvar)
  p2 <- which.min(pvar)

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
  } else if (method == "distcurva") {


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
