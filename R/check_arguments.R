.check_dats  <- function(dats)
{
  if ( sum(colnames(dats) %in% c("endpoint", "chemical", "concs")) != 3 )
  {
    rlang::abort("dataset does not have endpoint, chemical, and concs columns")
  } else
  {
    return(dats)
  }
}

.check_directionality <- function(directionality)
{
  if (sum(directionality %in% c(1, 0, -1)) != 1)
  {
    rlang::abort("only 1, 0, or -1 is allowed")
  } else
  {
    return(directionality)
  }
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
