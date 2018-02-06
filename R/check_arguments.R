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
      rlang::abort("threshold is a list but the names are not c(-1, 1)")
    }
  } else if (!is.numeric(as.numeric(threshold)))
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
