% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/linking_functions.R
\name{runLinking}
\alias{runLinking}
\title{Run Scale Linking}
\usage{
runLinking(data, method, verbose = FALSE, ...)
}
\arguments{
\item{data}{a \code{\linkS4class{PROsetta_data}} object. See \code{\link{loadData}} for loading a dataset.}

\item{method}{the linking method to use. Accepts:
\itemize{
  \item{\code{MM} for mean-mean method}
  \item{\code{MS} for mean-sigma method}
  \item{\code{HB} for Haebara method}
  \item{\code{SL} for Stocking-Lord method}
  \item{\code{FIXEDPAR} for fixed parameter calibration}
  \item{\code{CP} for calibrated projection using fixed parameter calibration on the anchor dimension}
  \item{\code{CPLA} for linear approximation of calibrated projection. This is identical to 'CP' but uses approximation in \code{\link{runRSSS}}}
  \item{\code{CPFIXEDDIM} for calibrated projection using mean and variance constraints on the anchor dimension}
}
Linear transformation methods (i.e., MM, MS, HB, SL) are performed with \code{\link[plink:plink-methods]{plink}} in \href{https://CRAN.R-project.org/package=plink}{'plink'} package.}

\item{verbose}{if \code{TRUE}, print status messages. (default = \code{FALSE})}

\item{...}{additional arguments to pass onto \code{\link[mirt]{mirt}} in \href{https://CRAN.R-project.org/package=mirt}{'mirt'} package.}
}
\value{
\code{\link{runLinking}} returns a \code{\link{list}} containing the scale linking results.
\itemize{
  \item{\code{constants}} linear transformation constants. Only available when linear transformation methods were used (i.e., MM, MS, HB, SL).
  \item{\code{ipar_linked}} item parameters calibrated to the response data, and linked to the metric of anchor item parameters.
  \item{\code{ipar_anchor}} anchor item parameters used in linking.
}
}
\description{
\code{\link{runLinking}} is a function for obtaining item parameters from the response data in the metric of supplied anchor item parameters.
}
\examples{
\donttest{
out_link <- runLinking(data_asq, "SL", technical = list(NCYCLES = 1000))
out_link$constants   # transformation constants
out_link$ipar_linked # item parameters linked to anchor
out_link <- runLinking(data_asq, "FIXEDPAR")
out_link$ipar_linked # item parameters linked to anchor
}
}
