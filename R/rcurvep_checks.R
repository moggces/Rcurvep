.check_dat_base  <- function(dat)
{
  cols_ori <- colnames(dat)
  cols_stand <- c("endpoint", "chemical", "conc", "resp")

  # remove extra columns
  cols_inter <- intersect(cols_ori, c(cols_stand, 'mask'))
  dat <- dat[, cols_inter]
  cols_removed <- setdiff(cols_ori, cols_inter)
  if (length(cols_removed) != 0) {
    rlang::warn(stringr::str_glue("{cols_removed} are removed in the input for run_rcurvep()"))
  }
  cols <- colnames(dat)

  # check columns
  if (!(all(rlang::has_name(dat, c(cols_stand, "mask"))) && length(cols) == 5))
  {
    if (!(all(rlang::has_name(dat, c(cols_stand))) && length(cols) == 4)) {
      rlang::abort("dataset needs to have endpoint, chemical, conc, resp, mask (optional) columns")
    }
  }

  dat2 <- dplyr::count(dat, .data$endpoint, .data$chemical, .data$conc)
  if (sum(dat2$n) != nrow(dat2))
  {
    rlang::abort("one concentration can have only one response")
  }

  endpoints <- unique(dat$endpoint)
  if (length(endpoints) > 1) {
    rlang::warn("same curvep parameters are applied to many endpoints")
  }

  return(dat)
}

.check_class <- function(obj, class_name, error_message) {
  if (class(obj) != class_name)
  {
    rlang::abort(error_message, class_name)
  }
  return(obj)
}

.check_mask_input <- function(vec, d) {
  if (is.null(vec)) return(vec)

  if (!all(vec == floor(vec))) {
    rlang::abort("input value must be integers")
  }
  min_n_conc <- d %>%
    dplyr::count(.data$endpoint, .data$chemical) %>%
    dplyr::pull(.data$n) %>%
    min(.)

  if (any(vec > min_n_conc)) {
    rlang::abort("input value must be smaller than total number of concs")
  }
  return(vec)
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

.check_dat  <- function(dat)
{
  cols <- colnames(dat)
  dev_cols <- c("endpoint", "chemical", "conc", "n_in", "N")
  beh_cols <- c("endpoint", "chemical", "conc", "resp")

  if (!(all(rlang::has_name(dat, c(beh_cols, "mask"))) && length(cols) == 5) &&
      !(all(rlang::has_name(dat, c(dev_cols, "mask"))) && length(cols) == 6) ) {
    if (!(all(rlang::has_name(dat, c(beh_cols))) && length(cols) == 4) &&
        !(all(rlang::has_name(dat, c(dev_cols))) && length(cols) == 5)) {
      rlang::abort("data is not a continuous/dichotomous dataset")
    }
  }
  return(dat)
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

.check_n_samples <- function(n_samples)
{
  if (!is.null(n_samples))
  {
    if (rlang::is_integer(as.integer(n_samples)) == FALSE || n_samples < 0)
    {
      rlang::abort("n_samples is not a valid number or is not NULL")
    }
  }
  return(n_samples)
}

.check_vdata <- function(vdata, dataset_type) {
  if (!all(!is.na(as.numeric(vdata)))) {
    rlang::abort("vdata is not NULL or all numeric")
  } else if (!is.null(vdata) && dataset_type == "dichotomous") {
    rlang::abort("currently numeric vdata is not supported for dichotomous dataset")
  }
  return(vdata)
}
