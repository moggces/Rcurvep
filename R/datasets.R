#' Zebrafish development endpoint data
#'
#' Data are processed from Biobide development plate incidence data
#'
#' @format a tibble with variables
#' \describe{
#' \item{endpoint}{toxicity endpoint}
#' \item{chemical}{chemical information}
#' \item{concs}{concentrations in log10(M) format}
#' \item{n_in}{number of incidence}
#' \item{N}{number of animals}
#' }
#'
#' @source Biobide study S-BBD-00016/15

"zfishdev"

#' Zebrafish behavior endpoint data
#'
#' Data are processed from Biobide behavior plate  data
#'
#' @format a tibble with variables
#' \describe{
#' \item{endpoint}{toxicity endpoint}
#' \item{chemical}{chemical information}
#' \item{concs}{concentrations in log10(M) format}
#' \item{resps}{normalized responses with baseline = 0}
#' }
#'
#' @source Biobide study S‐BBD‐0017/15

"zfishbeh"
