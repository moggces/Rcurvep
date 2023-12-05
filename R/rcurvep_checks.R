#_onerun
#_onerun(hill)
.check_dat_base  <- function(dat, check_one_concresp = TRUE)
{
  cols_ori <- colnames(dat)
  cols_stand <- c("endpoint", "chemical", "conc", "resp")

  # remove extra columns
  cols_inter <- intersect(cols_ori, c(cols_stand, 'mask'))
  dat <- dat[, cols_inter]
  cols_removed <- setdiff(cols_ori, cols_inter)
  if (length(cols_removed) != 0) {
    rlang::warn(stringr::str_glue("{cols_removed} are removed in the input"))
  }
  cols <- colnames(dat)

  # check columns
  if (!(all(rlang::has_name(dat, c(cols_stand, "mask"))) && length(cols) == 5))
  {
    if (!(all(rlang::has_name(dat, c(cols_stand))) && length(cols) == 4)) {
      rlang::abort(
        "dataset needs to have endpoint, chemical, conc, resp, mask (optional) columns"
      )
    }
  }

  if (check_one_concresp) {
    dat2 <- dplyr::count(dat, .data$endpoint, .data$chemical, .data$conc)
    if (sum(dat2$n) != nrow(dat2))
    {
      rlang::abort("one endpoint-chemical-concentration pair can have only one response. Use combi_run_rcurvep() instead")
    }
  }

  return(dat)
}

#_onerun
.check_mask_input <- function(vec, d) {

  if (rlang::has_name(d, "mask")) {
    rlang::warn("mask column exists; use original mask column")
    return(NULL)
  }
  if (any(vec == 0)) return(0)
  if (any(is.null(vec))) return(NULL)

  if (!all(vec == floor(vec))) {
    rlang::abort("input value must be integers")
  }
  if (any(vec < 0)) {
    rlang::abort("input value must be equal or larger than 0")
  }

  # calculate min of total number of concentrations
  min_n_conc <- d %>%
    dplyr::count(.data$endpoint, .data$chemical) %>%
    dplyr::pull(.data$n) %>%
    min()


  if (any(vec > min_n_conc)) {
    rlang::abort("input value must be smaller than total number of concs")
  }
  return(vec)
}

#_onerun
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


#_onerun
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

#_combirun
.check_config_name2 <- function(config = curvep_defaults(),  ...) {
  if (length(list(...)) == 0) {
    rlang::abort("Input parameters are needed, for example, TRSH = 5")
  }
  config <- .check_config_name(config = curvep_defaults(),  ...)
  return(config)
}


#_simulate, _combirun
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


#_simulate, _combirun
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

#_simulate, _combirun
.check_vdata <- function(vdata, dataset_type) {
  if (!all(!is.na(as.numeric(vdata)))) {
    rlang::abort("vdata is not NULL or all numeric")
  } else if (!is.null(vdata) && dataset_type == "dichotomous") {
    rlang::abort("currently numeric vdata is not supported for dichotomous dataset")
  }
  return(vdata)
}

#_onerun, _combirun
# .check_keep_sets <- function(keep_sets) {
#   keep_sets <- unique(keep_sets)
#   if (length(keep_sets) > 3 || !all(keep_sets %in% c("act_set", "resp_set", "fp_set"))) {
#     rlang::abort("Only a combination of act_set, resp_set, fp_set is allowed")
#   }
#   if (!'act_set' %in% keep_sets) {
#     rlang::abort("act_set is needed")
#   }
#   return(keep_sets)
# }

#_onerun, _combirun, _onerun(hill)
.check_keep_sets <- function(keep_sets, allowed_sets, must_set) {
  keep_sets <- unique(keep_sets)
  if (length(keep_sets) > length(allowed_sets) || !all(keep_sets %in% allowed_sets)) {
    sets <- stringr::str_c(allowed_sets, collapse = " ")
    rlang::abort(
      stringr::str_glue(
        "Only a combination {sets} is allowed.")
    )
  }
  if (!must_set %in% keep_sets) {
    rlang::abort(
      stringr::str_glue("{must_set} is needed.")
    )
  }
  return(keep_sets)
}

#_summary
.check_result_sets <- function(lsets) {
  if (!rlang::is_list(lsets) || length(lsets) > 3 || !all(names(lsets) %in% c("act_set", "resp_set", "fp_set"))) {
    rlang::abort("Only a list with names of act_set, resp_set, fp_set is allowed")
  }
  if (!"act_set" %in% names(lsets)) {
    rlang::abort("At least act_set tibble is needed")
  }
  return(lsets)
}

#_summary
.check_combirun_out <- function(combi_run_out) {
  if (!rlang::is_list(combi_run_out) || !rlang::has_name(combi_run_out, "result") || !rlang::has_name(combi_run_out, "config")) {
    rlang::abort("the combi_run_out is not a list or does not have result and config item.\n The input does not come from combi_run_rcurvep()")
  }
  return(combi_run_out)
}

#_bmr
.check_bmr_input <- function(d) {
  d <- .check_class(d, "rcurvep", "not a rcurvep object")
  if (!rlang::has_name(d$result$act_set, "sample_id")) {
    rlang::abort("The input act_set needs to have multiple samples (sample_id).
                 Please set n_samples in combi_run_rcurvep()")
  }
  return(d)
}

#_plot
.check_bmr_statsd <- function(d) {

  # adjust the colnames to the proper ones
  result <- d
  if (!rlang::has_name(d, "endpoint")) {
    col_names <- colnames(d)
    col_names <- c("TRSH", "pvar", col_names[-c(1,2)])
    result <- d %>% rlang::set_names(col_names)
    result <- result %>%
      dplyr::mutate(
        endpoint = "noname"
      ) %>%
      dplyr::select(
        .data$endpoint, dplyr::everything()
      )
  }
  return(result)
}

#_print, _plot, _bmr
.check_class <- function(obj, class_name, error_message) {
  if ( !class_name %in% class(obj) )
  {
    rlang::abort(error_message, class_name)
  }
  return(obj)
}

#_hillbase
.check_mask_on_concresp <- function(conc, resp, mask) {

  result <- list(Conc = conc, Resp = resp)

  len_inp <- purrr::map_dbl(list(conc, resp), length)
  if (length(unique(len_inp)) != 1 ) {
    rlang::abort("the length of Conc Resp is not the same.")
  }

  if (!is.null(mask)) {
    if (length(mask) != length(conc)) {
      rlang::abort("the length of Conc Resp Mask is not the same.")
    }
    result <- mask_resp_conc(conc, resp, mask)
    conc <- result$Conc
    resp <- result$Resp
  }

  if (length(conc) < 4) {
    rlang::abort("at least four resps are needed.")
  }

  return(result)
}

#_hillbase
# for hill model, f, theta, ui, ci can be modified
# for cc2 model, currently no additional parameters
.check_modl_args_exist <- function(ori, new) {
  if (!all(names(new) %in% names(ori))) {
    rlang::warn("The supplied arguments are not registered. No changes will be applied.")
  }
  config <- modifyList(ori, new)
  return(config)
}

#_hillbase
.check_modls_args_prefix <- function(args, modls) {

  syn <- paste0("(", stringr::str_c(modls, collapse = "|"), ")")
  if (!all(stringr::str_detect(names(args), syn))) {
    rlang::warn("Some supplied arguments do not have model type as the prefix. No changes will be applied.")
  }
  return(args)
}


#_hillbase
.check_modl_avail <- function(allowed, modls) {
  if (!all(modls %in% allowed)) {
    rlang::warn(
      stringr::str_glue(
        "{not_avail} not in the available model list.",
        not_avail = stringr::str_c(setdiff(modls, allowed), collapse = " ")
      )
    )
    result <- intersect(modls, allowed)
  } else {
    result <- modls
  }
  return(result)
}

#_hillbase
.check_modls_combi <- function(modls) {
  if (length(modls) > 1) {
    if ('cc2' %in% modls) {
      rlang::abort(">1 models are selected including cc2. Currently cc2 needs to be used alone.")
    }
  }
  return(modls)
}

#_run_fit
.check_modls_simu <- function(modls) {
  if (length(modls) > 1) {
    rlang::abort("Only hill or cc2 is allowed")
  } else if (!modls %in% c('hill', 'cc2')) {
    rlang::abort("Only hill or cc2 is allowed")
  }
  return(modls)
}

#_mergeobj
.check_objs_type <- function(objs) {

  # check the class
  objclass <- purrr::map_lgl(
    objs, ~ any(stringr::str_detect(class(.x), "rcurvep")))

  if (!all(objclass)) {
    rlang::abort("all objects need to be rcurvep objects")
  } else {

    # check the calculation
    objtype <- purrr::map_lgl(
      objs, ~ "config" %in% names(.x)
    )
    objsize <- purrr::map_int(
      objs, ~ length(.x$result)
    )

    if (!all(objtype)) {
      rlang::abort(
        "only results by curvep calculation are supported")
    }
    if (length(unique(objsize)) != 1) {
      rlang::abort(
        "different list size in results")
    }
  }

  return(objs)
}

#_summary
.check_inactivate <- function(inactivate) {
  inactivate <- na.omit(inactivate)
  if (!is.character(inactivate) & !rlang::is_integerish(inactivate) & !is.null(inactivate)) {
    rlang::abort("only string or integer is allowed for inactivate")
  }
  return(inactivate)
}
