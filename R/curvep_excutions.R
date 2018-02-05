create_curvep_input <- function(concs, resps, directionality, thres, mask, other_paras)
{
  rng <- 1000000
  if (directionality == -1) { rng <- -1000000 }

  other_paras['TRSH'] <- thres
  other_paras['RNGE'] <- rng

  if (sum(is.na(mask)) > 0)
  {
    #tibble here will be very slow
    result <- list(concs = list(concs), resps = list(resps), paras = list(other_paras))

  } else
  {
    result <- list(concs = list(concs), resps = list(resps), mask = list(mask), paras = list(other_paras))
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
  vals <- outp[!names(outp) %in% c('resp', 'corr', 'levels', 'Settings')]
  result <- vals %>% tibble::as_tibble(validate = FALSE)
  result[result == '-999'] <- NA
  return(result)
}


