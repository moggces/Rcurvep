% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tcplfit_hillbase.R
\name{fit_modls}
\alias{fit_modls}
\title{Fit one set of concentration-response data using types of models}
\usage{
fit_modls(Conc, Resp, Mask = NULL, modls = c("hill", "cnst"), ...)
}
\arguments{
\item{Conc}{A vector of log10 concentrations.}

\item{Resp}{A vector of numeric responses.}

\item{Mask}{Default = NULL or a vector of 1 or 0.
1 is for masking the respective response.}

\item{modls}{The model types for the fitting. Multiple values are allowed.
Currently Hill model (hill) and constant model (cnst) are implemented.
Default = c("hill", "cnst").}

\item{...}{The named input configurations for replacing the default configurations.
The input configuration needs to add model type as the prefix.
For example, hill_pdir = -1 will set the Hill fit only to the decreasing direction.}
}
\value{
A list of components named by the models.
The models are sorted by their AIC values. Thus, the first component has the best fit.
\subsection{hill}{

Fit output from Hill equation
\itemize{
\item modl: model type, i.e., hill
\item fit: fittable, 1 (yes) or 0 (no)
\item aic: AIC value
\item tp: model top, <0 means the fit for decreasing direction is preferred
\item ga: ac50 (log10 scale)
\item gw: Hill coefficient
\item er: scale term for Student's t distribution
}
}

\subsection{cnst}{

Fit output from constant model
\itemize{
\item modl: model type, i.e., cnst
\item fit: fittable?, 1 or 0
\item aic: AIC value
\item er: scale term
}
}
}
\description{
A convenient function to fit data using available models
and to sort the outcomes by AIC values.
}
\details{
The backbone of fit using hill and cnst is based on the implementation from tcpl package.
But the lower bound of ga is lower by log10(1/100).
}
\examples{

concd <- c(-9, -8, -7, -6, -5, -4)
respd <- c(0, 2, 30, 40, 50, 60)
maskd <- c(0, 0, 0, 0, 0, 1)

# run hill only
fit_modls(concd, respd, modls = "hill")

# run hill only + increasing direction only
fit_modls(concd, respd, modls = "hill", hill_pdir = 1)

# run with mask at the highest concentration
fit_modls(concd, respd, maskd)


}
\seealso{
\code{\link[tcpl:Models]{tcpl::tcplObjHill()}}, \code{\link[tcpl:Models]{tcpl::tcplObjCnst()}}, \code{\link[=get_hill_fit_config]{get_hill_fit_config()}}
}