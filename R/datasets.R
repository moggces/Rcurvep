#' Subsets of concentration response datasets from zebrafish developmental toxicity assays
#'
#' The datasets contain 4 toxicity endpoints and 3 chemicals.
#'
#'
#' @format A tibble with 96 rows and 5 columns:
#'
#' \describe{
#'   \item{endpoint}{endpoint name + at time point measured}
#'   \item{chemical}{chemical name + CASRN}
#'   \item{conc}{concentrations in log10(M) format}
#'   \item{n_in}{number of incidence}
#'   \item{N}{number of embryos}
#' }
#'
#' @source Biobide study S-BBD-00016/15

"zfishdev"

#' Subsets of concentration response datasets from zebrafish neurotoxicity assays
#'
#' The datasets contain 11 toxicity endpoints and 2 chemicals.
#' The responses have been normalized so that the baseline is 0.
#'
#' @format A tibble with 2123 rows and 4 columns:
#' \describe{
#'   \item{endpoint}{endpoint name}
#'   \item{chemical}{chemical name + CASRN}
#'   \item{conc}{concentrations in log10(M) format}
#'   \item{resp}{responses after normalized using the vehicle control on each plate}
#' }
#'
#' @source Biobide study S-BBD-0017/15

"zfishbeh"

#' Full sets of concentration response datasets from zebrafish developmental toxicity assays
#'
#' The datasets contain 4 toxicity endpoints and 32 chemicals.
#'
#'
#' @format A tibble with 512 rows and 5 columns:
#' @seealso [zfishdev]
#'
#' @source Biobide study S-BBD-00016/15

"zfishdev_all"


#' Activity output based on simulated datasets using zfishdev_all dataset
#'
#' The data is an rcurvep object from the [combi_run_rcurvep()].
#' See [combi_run_rcurvep()] for the code to reproduce this dataset.
#'
#' @format A list of two named components: result and config.
#'   The result component is a list with one component: act_set.
#'
#' @seealso [estimate_dataset_bmr()]
#'

"zfishdev_act"



