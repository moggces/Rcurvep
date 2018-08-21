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
identify_basenoise_threshold <- function(df, endpoint, chemical, threshold, direction, potency, p1 = NULL, p2 = NULL, plot = TRUE, n_endpoint_page = 4) {


  #id <- rlang::syms(id)
  endpoint <- rlang::sym(endpoint)
  direction <- rlang::sym(direction)
  chemical <- rlang::sym(chemical)
  threshold <- rlang::sym(threshold)
  potency <- rlang::sym(potency)

  thres_vars_pre <- df %>%
    dplyr::group_by(!!endpoint, !!direction, !!chemical, !!threshold) %>%
    dplyr::summarize(sd_pod = sd(!!potency), n_rep = n()) %>%
    dplyr::ungroup() %>%
    dplyr::group_by(!!endpoint, !!direction, !!threshold) %>%
    dplyr::summarize(pooled_variance  = sum((.data$sd_pod^2)*(.data$n_rep - 1))/(sum(.data$n_rep) - n())) %>%
    dplyr::ungroup()


   thres_vars <- thres_vars_pre %>%
     tidyr::nest(-c(!!endpoint, !!direction)) %>%
     dplyr::mutate(
       outd = purrr::map(data, function(x) cal_dist_curva(x[[as.character(threshold)]], x$pooled_variance, p1 = p1, p2 = p2))
     ) %>%
     tidyr::unnest()

   select_thres <- thres_vars %>%
    tidyr::nest(-c(!!endpoint, !!direction)) %>%
    dplyr::mutate(
      thresCurva = purrr::map(data, function(x) pick_threshold(x[[as.character(threshold)]], x$dist2l, x$curvature, method = "distcurva")),
      thresDist = purrr::map(data, function(x) pick_threshold(x[[as.character(threshold)]], x$dist2l, method = "dist"))

    ) %>%
    tidyr::unnest()

   if (plot == TRUE)
   {
     select_thres_g <- select_thres %>%
       dplyr::select(!!endpoint, !!direction, !!threshold, "pooled_variance", "dist2l", "curvature") %>%
       tidyr::gather(type, value, c("pooled_variance", "dist2l", "curvature")) %>%
       dplyr::mutate(type = ordered(type, levels =  c("pooled_variance", "dist2l", "curvature") )) %>%
       dplyr::mutate(endpoint = !!endpoint, direction = !!direction, threshold = !!threshold) %>%
       dplyr::mutate(direction = as.factor(direction))

     d <- select_thres_g %>% dplyr::pull(endpoint) %>% unique
     if (length(d) < n_endpoint_page) n_endpoint_page <- length(d)
     d <- split(d, ceiling(seq_along(d)/n_endpoint_page))

     pl <- purrr::map(d, function(x) {
       #outthres_gp <- select_thres_g %>% dplyr::filter(!!rlang::expr((!!endpoint) %in% x))
       outthres_gp <- select_thres_g %>% dplyr::filter(endpoint %in% x)
       p <- ggplot2::ggplot(outthres_gp,
                  #ggplot2::aes_(x = threshold, y = "value", color = as.factor(direction)))
                  ggplot2::aes(x = threshold, y = value, color = direction))
       p <- p + ggplot2::geom_point(size = 2) +
         ggplot2::geom_line() +
         ggplot2::facet_grid(type ~ endpoint, scales = "free")
       return(p)
     })
     print(pl)
   }

  return(select_thres)
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
