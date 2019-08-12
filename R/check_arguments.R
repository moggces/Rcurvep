
.check_dat_base  <- function(dats)
{
  cols <- colnames(dats)
  cols_stand <- c("endpoint", "chemical", "conc", "resp")
  if (length(setdiff(cols, cols_stand)) != 0)
  {
    rlang::abort("dataset needs to have endpoint, chemical, conc, resp columns")
  }

  dats2 <- dplyr::count(dats, endpoint, chemical, conc)
  if (sum(dats2$n) != nrow(dats2))
  {
    rlang::abort("one concentration can have only one response")
  }

  endpoints <- unique(dats$endpoint)
  if (length(endpoints) > 1) {
    rlang::warn("same curvep parameters are applied to many endpoints")
  }

  return(dats)
}

.check_directionality <- function(directionality)
{
  if (length(directionality) != 1 | sum(directionality %in% c(1, -1)) != 1)
  {
    rlang::abort("only 1, -1 is allowed")
  } else
  {
    return(directionality)
  }
}

.check_class <- function(obj, class_name, error_message) {
  if (class(obj) != class_name)
  {
    rlang::abort(error_message, class_name)
  }
  return(obj)
}

.check_config_name <- function(config = curvep_defaults(),  ...) {

  config <- .check_class(config, 'curvep_config', "config is absent or corrupt")

  args <- list(...)
  config <- modifyList(config, args)

  defaults <- curvep_defaults()
  para_diff <- setdiff(names(config), names(defaults))

  if (length(para_diff) > 0) {
    rlang::abort("nonexistent curvep parameters are added")
  }

  return(config)
}

.check_config_value <- function(config) {


  if (rlang::is_double(config$TRSH) == FALSE | config$TRSH < 0)
  {
    rlang::abort("TRSH has to be a positive double value")
  }
  if (rlang::is_double(config$RNGE) == FALSE)
  {
    rlang::abort("RNGE has to be a double value")
  }
  if (rlang::is_double(config$MXDV) == FALSE | config$MXDV < 0 ) {
    rlang::abort("MXDV has to be a positive double value")
  }
  if (rlang::is_integer(config$BSFT) | config$BSFT < 0) {
    rlang::abort("BSFT has to be a positive integer value")
  }
  if (rlang::is_integer(config$USHP) | config$USHP < 0) {
    rlang::abort("USHP has to be a positive integer value")
  }
  if (rlang::is_logical(config$TrustHi) == FALSE | rlang::is_logical(config$StrictImp) == FALSE) {
    rlang::abort("TrustHi and TrustHi have to be a TRUE/FALSE value")
  }

  return(config)



}

.check_n_sample <- function(n_sample)
{
  if (!is.null(n_sample))
  {
    if (rlang::is_integer(rlang::as_integer(n_sample)) == FALSE)
    {
      rlang::abort("n_sample is not a valid number or is not NULL")
    }
  }
  return(n_sample)

}



.check_dats  <- function(dats)
{
  cols <- colnames(dats)
  dev_cols <- c("endpoint", "chemical", "conc", "n_in", "N")
  beh_cols <- c("endpoint", "chemical", "conc", "resp")

  if (length(setdiff(cols, dev_cols)) != 0 | length(setdiff(cols, beh_cols)) != 0)
  {
    rlang::abort("data is not a continuous/dichotomous dataset")
  } else
  {
    return(dats)
  }
}

.check_dats_reparam  <- function(dats)
{
  if (
    sum(colnames(dats) %in% c("threshold", "endpoint", "chemical", "direction", "repeat_id", "input")) != 6 )
  {
    rlang::abort("dataset is not the unsimplified output from run_curvep_job()")
  } else
  {
    return(dats)
  }
}



.check_directionality_reparam <- function(directionality)
{
  if (!is.null(directionality)) {
    if (length(directionality) != 1 | sum(directionality %in% c(1,-1)) != 1) {
      rlang::abort("only 1,-1, NULL is allowed")
    }
  }
  return(directionality)
}



.check_threshold <- function(threshold, directionality)
{
  directionality <- .check_directionality(directionality)
  if (rlang::is_list(threshold))
  {
    if ( sum(names(threshold) %in% c(1, -1)) != length(threshold) )
    {
      rlang::abort("threshold is a list but the names are not 1 or -1")
    }
    if (directionality == 0)
    {
      if (sum(names(threshold) %in% c(1, -1)) != 2)
      {
        rlang::abort("directionality = 0 but the named list threshold does not have 1 and -1")
      }
    } else if (directionality == 1)
    {
      if (sum(names(threshold) %in% c(1)) != 1)
      {
        rlang::abort("directionality = 1 but the named list threshold does not have 1")
      }
    } else if (directionality == -1)
    {
      if (sum(names(threshold) %in% c(-1)) != 1)
      {
        rlang::abort("directionality = -1 but the named list threshold does not have -1")
      }
    }
    if (sum(purrr::map_lgl(threshold, is.numeric)) != length(threshold) )
    {
      rlang::abort("threshold is not a numeric vector")
    }
  } else if (!is.numeric(threshold))
  {
    rlang::abort("threshold is not a numeric vector")
  }
  return(threshold)
}

.check_threshold_reparam <- function(threshold) {

  if (!is.null(threshold)) {
    if (is.list(threshold) | length(threshold) != 1 | !is.numeric(threshold) | threshold < 0) {
      rlang::abort("only allow one numeric positive threshold")
    }
  }
  return(threshold)
}

.check_other_paras <- function(other_paras)
{
  allowed_paras <- c('MXDV', 'CARR', 'BSFT', 'USHP', 'TrustHi', 'StrictImp')
  if (!rlang::is_list(other_paras))
  {
    rlang::abort("other_paras is not a list")
    if (rlang::is_names(other_paras))
    {
      if (sum(names(other_paras) %in% allowed_paras) != length(other_paras))
      {
        rlang::abort("other_paras contains unknown parameters. only `allowed_paras`")
      }
    }
  }
  return(other_paras)
}

.check_plot_firstin <- function(x) {
  if (inherits(x, "rcurvep_thres_stats"))
  {
    df <- x[["stats"]]

  } else {
    df <- x

  }
}

.check_plot_paras <- function(...) {
  dots <- list(...)

  endpoint <- ifelse(!is.null(dots$endpoint), dots$endpoint, "endpoint")
  direction <- ifelse(!is.null(dots$direction), dots$direction, "direction")
  threshold <- ifelse(!is.null(dots$threshold), dots$threshold, "threshold")
  n_endpoint_page <- ifelse(!is.null(dots$n_endpoint_page), dots$n_endpoint_page, 4)

  return(list(df = df, endpoint = endpoint, direction = direction, threshold = threshold, n_endpoint_page = n_endpoint_page))
}
