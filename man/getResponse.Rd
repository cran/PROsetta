% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions.R
\name{getResponse}
\alias{getResponse}
\title{Extract scale-wise response}
\usage{
getResponse(d, scale_id = "all", person_id = FALSE)
}
\arguments{
\item{d}{a \code{\linkS4class{PROsetta_data}} object.}

\item{scale_id}{scale IDs to extract response. If \code{all}, use all scale IDs. (default = \code{all})}

\item{person_id}{if \code{TRUE}, also return person IDs. (default = \code{FALSE})}
}
\value{
\code{\link{getResponse}} returns a \code{\link{data.frame}} containing scale-wise response.
}
\description{
\code{\link{getResponse}} is a helper function for extracting scale-wise response from a \code{\linkS4class{PROsetta_data}} object.
}
\examples{
getResponse(data_asq)
getResponse(data_asq, 1)
getResponse(data_asq, 2)
getResponse(data_asq, c(1, 2))
getResponse(data_asq, c(2, 1))
getResponse(data_asq, c(1, 2), person_id = TRUE)

}
