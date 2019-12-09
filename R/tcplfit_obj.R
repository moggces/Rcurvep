
ObjHillnorm <-  function(p, lconc, resp) {

  mu <- p[1]/(1 + 10^((p[2] - lconc)*p[3]))
  sum(dnorm(x=resp,mean=mu,sd=p[4],log = TRUE))

}

tcplHillVal <- function(logc, tp, ga, gw, bt = 0) {

  bt + (tp - bt)/(1 + 10^((ga - logc)*gw))

}

tcplHillConc <- function(val, tp, ga, gw, bt = 0) {

  ga - log10((tp - bt)/(val - bt) - 1)/gw

}

tcplHillACXX <- function(XX, tp, ga, gw, bt = 0) {

  y <- tp * XX/100
  ga - log10((tp - bt)/(y - bt) - 1)/gw

}

#-------------------------------------------------------------------------------
# tcplObjCnst: Generate a constant model objective function to optimize
#-------------------------------------------------------------------------------

#' @rdname Models
#'
#' @section Constant Model (cnst):
#' \code{tcplObjCnst} calculates the likelyhood for a constant model at 0. The
#' only parameter passed to \code{tcplObjCnst} by \code{p} is the scale term
#' \eqn{\sigma}. The constant model value \eqn{\mu_{i}}{\mu[i]} for the
#' \eqn{i^{th}}{ith} observation is given by:
#' \deqn{\mu_{i} = 0}{\mu[i] = 0}
#'
#' @importFrom stats dt
#' @export

tcplObjCnst <- function(p, resp) {

  ### This function takes creates an objective function to be optimized using
  ### the starting constant model parameter, and response.
  ###
  ### Arguments:
  ###   p:     a numeric vector of length 1 containg the starting values for
  ###          the constant model, in order: log error term
  ###   lresp: a numeric vector containing the response values to produce the
  ###          objective function
  ###
  ### Value:
  ###   An objective function for the constant model and the given resp data

  mu <- 0
  sum(dt((resp - mu)/exp(p[1]), df = 4, log = TRUE) - p[1])

}

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# tcplObjHill: Generate a hill model objective function to optimize
#-------------------------------------------------------------------------------

#' @rdname Models
#'
#' @section Hill Model (hill):
#' \code{tcplObjHill} calculates the likelyhood for a 3 parameter Hill model
#' with the bottom equal to 0. The parameters passed to \code{tcplObjHill} by
#' \code{p} are (in order) top (\eqn{\mathit{tp}}), log AC50 (\eqn{\mathit{ga}}), hill
#' coefficient (\eqn{\mathit{gw}}), and the scale term (\eqn{\sigma}). The hill model
#' value \eqn{\mu_{i}}{\mu[i]} for the \eqn{i^{th}}{ith} observation is given
#' by:
#' \deqn{
#' \mu_{i} = \frac{tp}{1 + 10^{(\mathit{ga} - x_{i})\mathit{gw}}}
#' }{
#' \mu[i] = tp/(1 + 10^(ga - x[i])*gw)}
#' where \eqn{x_{i}}{x[i]} is the log concentration for the \eqn{i^{th}}{ith}
#' observation.
#'
#' @importFrom stats dt
#' @export

tcplObjHill <- function(p, lconc, resp) {

  ### This function takes creates an objective function to be optimized using
  ### the starting hill parameters, log concentration, and response.
  ###
  ### Arguments:
  ###   p:     a numeric vector of length 4 containg the starting values for
  ###          the hill model, in order: top, log AC50, hill
  ###          coefficient, and log error term
  ###   lconc: a numeric vector containing the log concentration values to
  ###          produce the objective function
  ###   lresp: a numeric vector containing the response values to produce the
  ###          objective function
  ###
  ### Value:
  ###   An objective function for the hill model and the given conc-resp data

  mu <- p[1]/(1 + 10^((p[2] - lconc)*p[3]))
  sum(dt((resp - mu)/exp(p[4]), df = 4, log = TRUE) - p[4])

}

#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# tcplObjGnls: Generate a gain-loss model objective function to optimize
#-------------------------------------------------------------------------------

#' @rdname Models
#'
#' @section Gain-Loss Model (gnls):
#' \code{tcplObjGnls} calculates the likelyhood for a 5 parameter model as the
#' product of two Hill models with the same top and both bottoms equal to 0.
#' The parameters passed to \code{tcplObjGnls} by \code{p} are (in order) top
#' (\eqn{\mathit{tp}}), gain log AC50 (\eqn{\mathit{ga}}), gain hill coefficient (\eqn{gw}),
#' loss log AC50 \eqn{\mathit{la}}, loss hill coefficient \eqn{\mathit{lw}}, and the scale
#' term (\eqn{\sigma}). The gain-loss model value \eqn{\mu_{i}}{\mu[i]} for the
#' \eqn{i^{th}}{ith} observation is given by:
#' \deqn{
#' g_{i} = \frac{1}{1 + 10^{(\mathit{ga} - x_{i})\mathit{gw}}}
#' }{
#' g[i] = 1/(1 + 10^(ga - x[i])*gw)}
#' \deqn{
#' l_{i} = \frac{1}{1 + 10^{(x_{i} - \mathit{la})\mathit{lw}}}
#' }{
#' l[i] = 1/(1 + 10^(x[i] - la)*lw)}
#' \deqn{\mu_{i} = \mathit{tp}(g_{i})(l_{i})}{\mu[i] = tp*g[i]*l[i]}
#' where \eqn{x_{i}}{x[i]} is the log concentration for the \eqn{i^{th}}{ith}
#' observation.
#'
#' @importFrom stats dt
#' @export

tcplObjGnls <- function(p, lconc, resp) {

  ### This function takes creates an objective function to be optimized using
  ### the starting gain-loss parameters, log concentration, and response.
  ###
  ### Arguments:
  ###   p:     a numeric vector of length 5 containg the starting values for
  ###          the gain-loss model, in order: top, gain log AC50, gain hill
  ###          coefficient, loss log AC50, loss hill coefficient and log error
  ###          term
  ###   lconc: a numeric vector containing the log concentration values to
  ###          produce the objective function
  ###   lresp: a numeric vector containing the response values to produce the
  ###          objective function
  ###
  ### Value:
  ###   An objective function for the gain-loss model and the given conc-resp
  ###   data

  gn <- 1/(1 + 10^((p[2] - lconc)*p[3]))
  ls <- 1/(1 + 10^((lconc - p[4])*p[5]))
  mu <- p[1]*gn*ls
  sum(dt((resp - mu)/exp(p[6]), df = 4, log = TRUE) - p[6])

}
