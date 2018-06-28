
create_curvep_input <- function(d, threshold, paras) {

  rng <- 1000000
  if (unique(d$directionality_u) == -1) { rng <- -1000000 }

  paras['TRSH'] <- threshold
  paras['RNGE'] <- rng

  if (sum(is.na(d$mask)) > 0)
  {
    #tibble here will be very slow
    result <- list(concs = list(d$concs), resps = list(d$resps), paras = paras )

  } else
  {
    result <- list(concs = list(d$concs), resps = list(d$resps), mask = list(d$mask), paras = paras)
  }
  return(result)
}


fill_directions <- function(dd, directionality) {

  if (directionality == 0)
  {
    result <- lapply(list(-1, 1), function(dirn){
      dd %>%
        dplyr::mutate(directionality_u = dirn)
    }) %>% dplyr::bind_rows()
  } else
  {
    result <- dd
    result <- result %>% dplyr::mutate(directionality_u  = directionality)
  }

  return(result)
}


run_curvep <- function(x) {
  paras <- unlist(x$paras, recursive = FALSE)
  concs <- unlist(x$concs)
  resps <- unlist(x$resps)
  valid <- !is.na(concs)
  mask <- x[['mask']] #the null
  result <- do.call(curvep, c(list(concs[valid]), list(resps[valid]), mask, paras)) #call curvep , c() to combine list
  return(result)
}


tabulate_curvep_output <- function(x) {

  outp <- x
  vals <- outp[!names(outp) %in% c('resp', 'corr', 'xx', 'ECxx', 'Cxx','Settings')]
  result <- vals %>% tibble::as_tibble(validate = FALSE)
  result[result == '-999'] <- NA
  return(result)
}


