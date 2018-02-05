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

  result <- ribble::tibble(
    dist2l = dists,
    der2 = spl_der2$y,
    der1 = spl_der1$y,
    curvature = spl_der2$y/(((1 + (spl_der1$y)^2))^3/2),
    p1 = p1, p2 = p2)

  return(result)

}


suggest_BMR <- function(thres, dist2l, curvature = NULL, method = c("dist", "distcurva")) {
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
