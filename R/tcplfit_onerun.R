# need to report unknown parameters
# ... needs hill_pdir
run_fit <- function(d, modls = c("hill", "cnst"), keep_sets = c('fit_set', 'resp_set'), n_samples = NULL, ...) {

  args <- list(...)

  nestd <- nest_fit_dataset(d, nest_cols = c("endpoint", "chemical"))
  fitd <- cal_fit_dataset(nestd = nestd, modls = modls, args = args, fit_type = "ori")

  if(!is.null(n_samples)) {
    simud <- create_hillsimu_dataset(fitd = fitd, n_samples = n_samples, pdir = args$hill_pdir)
    nestd <- nest_fit_dataset(simud, nest_cols = c("endpoint", "chemical", "sample_id"))
    fitd <- cal_fit_dataset(nestd = nestd, modls = modls, args = args, fit_type = "hill_simu")
  }

  fitd_clean <- clean_fit_output(fitd = fitd, keep_sets = keep_sets)
  sets <- merge_fit_output(fitd_clean, keep_sets = keep_sets)

  result <- list(result = sets, result_nested = fitd)

  return(result)
}

merge_fit_output <- function(fitd, keep_sets) {

  tbl_names <- list(
    fit_set = c("out_paras"),
    resp_set = c("out_resps"),
    act_set = c("in_summary", "activity")
  )

  result <- purrr::map(
    keep_sets, function(x, tbl_names, obj)
      suppressWarnings(tidyr::unnest(
        obj %>% dplyr::select(tidyselect::one_of(c("endpoint", "chemical", "sample_id", tbl_names[[x]])))
      )),
    tbl_names = tbl_names,
    obj = fitd
  ) %>% rlang::set_names(keep_sets)
  return(result)

}

summarize_fit_output <- function(d, thr_resp = 20, perc_resp = 10, ci_level = 0.95, extract_only = FALSE) {

  fitd <- d$result_nested
  fitd_clean <- clean_fit_output(fitd, keep_sets = "act_set", thr_resp = thr_resp, perc_resp = perc_resp)
  act_set <- merge_fit_output(fitd_clean, keep_sets = c("act_set"))
  d[['result']] <- c(d$result, act_set)

  if (!extract_only) {
    lsets <- d$result
    base_cols <- get_base_cols(lsets)
    lsets_n <- get_nested_joined_sets(lsets, base_cols)
    lsets_n <- make_act_na_highconc(lsets_n)
    lsets_n <- make_act_na_0(lsets_n)
    lsets_n <- summarize_actsets(lsets_n, ci_level = ci_level)
    d[['act_summary']] <- unnest_joined_sets(lsets_n, base_cols, "act_summary")
  }

  return(d)
}

make_act_na_0 <- function(nestd) {
  result <- nestd %>%
    dplyr::mutate(
      act_set = purrr::map(.data$act_set, make_act_na_0_in)
    )
  return(result)
}

make_act_na_0_in <- function(act_set) {

  potency_cols <- c("slope", "Emax")
  suppressWarnings(result <- act_set %>%
                     dplyr::mutate_at(
                       dplyr::vars(tidyselect::one_of(potency_cols)),
                       ~ replace(.x, is.na(.x), 0)
                     ))

  return(result)
}

cal_fit_dataset <- function(nestd, modls, args, fit_type = c("ori", "hill_simu")) {
  if (fit_type == "ori") {
    result <- nestd %>%
      dplyr::mutate(
        output = purrr::map(
          .data$input, ~ do.call(fit_modls,c(list(Conc = .x$conc, Resp = .x$resp, Mask = NULL, modls = modls), args))
        )
      )
  } else if (fit_type == "hill_simu") {

    args[['hill_pdir']] <- NULL

    result <- nestd %>%
      dplyr::mutate(
        output = purrr::map(
          .data$input, ~ do.call(fit_modls,c(list(Conc = .x$conc, Resp = .x$resp, hill_pdir = unique(.x$direction), Mask = NULL, modls = modls), args))
        )
      )
  }

  return(result)
}

clean_fit_output <- function(fitd, keep_sets, thr_resp = NULL, perc_resp = NULL) {
  result <- suppressMessages(
    purrr::reduce(
      purrr::map(
        keep_sets, clean_fit_output_in,
        nestd = fitd, thr_resp = thr_resp, perc_resp = perc_resp),
      dplyr::left_join)
  )
  return(result)
}

clean_fit_output_in <- function(nestd, keep_set, thr_resp = NULL, perc_resp = NULL) {

  if (keep_set == "fit_set") {
    result <- nestd %>%
      dplyr::mutate(
        out_paras = purrr::map(.data$output, extract_fit_para)
      )
  }
  if (keep_set == 'resp_set') {
    result <- nestd %>%
      dplyr::mutate(
        out_resps = purrr::map2(.data$input, .data$output, extract_fit_resp)
      )
  }
  if (keep_set == "act_set") {
    result <- nestd %>%
      dplyr::mutate(
        in_summary = purrr::map(.data$input, extract_input_summary),
        activity = purrr::map2(.data$input, .data$output, extract_fit_activity, thr_resp = thr_resp, perc_resp = perc_resp)
      )
  }
  result <- result %>% dplyr::select(-.data$input, -.data$output)

  return(result)
}

nest_fit_dataset <- function(d, nest_cols) {

  nest_colsq <- rlang::syms(c(nest_cols, 'conc'))

  result <- d %>%
    dplyr::arrange(!!!nest_colsq) %>%
    tidyr::nest(-tidyselect::one_of(nest_cols), .key = "input")

  return(result)
}

cal_hillfit_dataset_in <- function(nestd, pdir, ...) {

  args <- list(...)

  if (!is.null(pdir)) {
    result <- nestd %>%
      dplyr::mutate(
        output = purrr::map(
          .data$input, ~ do.call(
              fit_modls,
              c(list(Conc = .x$conc, Resp = .x$resp, Mask = NULL, types = "hill", pdir = pdir), args)
          )
        )
      )
  } else { # direction column in the input
    result <- nestd %>%
      dplyr::mutate(
        output = purrr::map(
          .data$input,  ~ do.call(
              fit_modls,
              c(list(Conc = .x$conc, Resp = .x$resp, Mask = NULL, types = "hill", pdir = unique(.x$direction)), args)
            )
        )
      )
  }
  return(result)

}

create_hillsimu_dataset <- function(fitd, n_samples, pdir) {
  result <- fitd %>%
    dplyr::mutate(
      simud = purrr::map2(
        .data$input, .data$output,
        ~ dplyr::bind_rows(
          replicate(n_samples, create_hillsimu_resp(.x, .y), simplify = FALSE),
          .id = "sample_id"))
    )

  # create direction column and remove the original fit input and output
  result <- result %>%
    dplyr::mutate(
      direction = purrr::map_dbl(.data$output, select_direction, pdir = pdir)
    ) %>%
    dplyr::select(-.data$input, -.data$output) %>%
    tidyr::unnest()

  return(result)
}

cal_hillfit_dataset <- function(d, pdir,  n_samples, ...) {

  d <- nest_fit_dataset(d, nest_cols = c("endpoint", "chemical"))

  fitd <- cal_hillfit_dataset_in(nestd = d, pdir = pdir, ...)
  result <- fitd

  if (!is.null(n_samples)) {
    simud <- create_simu_dataset(fitd = fitd, pdir = pdir, n_samples = n_samples)
    simud <- nest_fit_dataset(simud, nest_cols = c("endpoint", "chemical", "sample_id"))
    result <- cal_hillfit_dataset_in(nestd = simud, pdir = NULL, ...)
  }

  return(result)
}

select_direction <- function(out, pdir) {

  tp <- out[['hill']]$tp

  if (!is.null(pdir)) {
    if (length(pdir) == 1) {
      return(pdir)
    }
  }

  result <- 1
  if(!is.na(tp)) {
    if (tp < 0) result <- -1
  }
  return(result)

}

tibble_fit_para <- function(modl) {
  pars <- tibble::as_tibble(modl)
  modl_type <- pars$modl
  result <- pars %>%
    dplyr::select(-.data$modl) %>%
    dplyr::rename_all(.funs = list(~ stringr::str_c(modl_type, "_", .)))
  return(result)
}

extract_fit_para <- function(lmodl) {
  result <- purrr::map_dfc(lmodl, tibble_fit_para)
  result[['win_modl']] <- lmodl[[1]]$modl
  return(result)
}

extract_fit_resp <- function(inp, out) {
  result <- inp
  win_modl <- out[[1]]
  fit_resps <- rep(as.numeric(NA), length(inp$conc))
  if (win_modl$modl == "cnst") {
    fit_resps <- rep(median(inp$resp), length(inp$conc))
  } else if (win_modl$modl == "hill") {
    fit_resps <- tcplHillVal(inp$conc, win_modl$tp, win_modl$ga, win_modl$gw)
  }
  result[['fitted_resp']] <- fit_resps
  return(result)
}


create_hillsimu_resp <- function(inp, out) {
  result <- inp
  yhat <- tcplHillVal(inp$conc, out[['hill']]$tp, out[['hill']]$ga, out[['hill']]$gw)
  resid <- inp$resp - yhat
  sampleresids <- sample(resid,length(resid),replace = TRUE)
  this_y <- yhat + sampleresids
  result$resp <- this_y
  return(result)
}






extract_fit_activity <- function(inp, out, thr_resp, perc_resp) {

  win_modl <- out[[1]]
  med_resp <- median(inp$resp)
  EC50 <- NA_real_
  Emax <- NA_real_
  POD <- NA_real_
  ECxx <- NA_real_
  slope <- NA_real_
  hit <- 0

  if (win_modl$modl == "cnst") {
    if (abs(med_resp) > thr_resp) hit <- 1

  } else if (win_modl$modl == "hill") {
    Emax <- win_modl$tp
    EC50 <- win_modl$ga
    slope <- win_modl$gw
    if (!is.na(Emax)) {
      ECxx <- tcplHillACXX(perc_resp, win_modl$tp, win_modl$ga, win_modl$gw)
      if (Emax < 0) thr_resp <- thr_resp*-1
      if (abs(Emax) > abs(thr_resp)) POD <- tcplHillConc(thr_resp, win_modl$tp, win_modl$ga, win_modl$gw)
      if (!is.na(POD) && POD < max(inp$conc)) hit <- 1
    }
  }

  result <- tibble::tibble(
    EC50 = EC50,
    slope = slope,
    Emax = Emax,
    ECxx = ECxx,
    POD = POD,
    hit = hit
    )
  return(result)

}


extract_hill_activity <- function(out, bmr) {
  if (out$tp < 0) bmr <- bmr*-1

  result <- list(
    EC50 = out[['ga']],
    Emax = out[['tp']],
    POD = suppressWarnings(tcplHillConc(bmr, out$tp, out$ga, out$gw))
  ) %>% tibble::as_tibble()

  return(result)
}

extract_hillfit_stats <- function(inp, out) {
  yfit <- tcplHillVal(inp$conc, out$tp, out$ga, out$gw)
  result <- tibble::tibble(
    r2 = cor(inp$resp, yfit)^2,
    rmse = sqrt(mean((yfit - inp$resp)^2, na.rm = TRUE))
  )
  return(result)
}

