% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/core_functions.R
\name{computeResponseProbability}
\alias{computeResponseProbability}
\title{(internal) compute response probability}
\usage{
computeResponseProbability(ipar, model, theta_grid)
}
\arguments{
\item{ipar}{a \code{\link{data.frame}} containing item parameters.}

\item{model}{the item model to use. Accepts \code{grm} or {gpcm}.}

\item{theta_grid}{theta values to compute probability values at.}
}
\value{
\code{\link{computeResponseProbability}} returns an item-wise list of probability matrices.
}
\description{
\code{\link{computeResponseProbability}} is an internal function for computing response probability from a set of item parameters.
}
\examples{
ipar <- PROsetta:::extractAnchorParameters(data_asq, FALSE)
theta_q <- seq(-4, 4, .1)
p <- PROsetta:::computeResponseProbability(ipar, "grm", theta_q)

plot(
  0, 0, type = "n", xlim = c(-4, 4), ylim = c(0, 1),
  xlab = "Theta", ylab = "Response probability"
)
lines(theta_q, p[[1]][, 1])
lines(theta_q, p[[1]][, 2])
lines(theta_q, p[[1]][, 3])
lines(theta_q, p[[1]][, 4])
lines(theta_q, p[[1]][, 5])

}
\keyword{internal}
