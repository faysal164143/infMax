% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/infMaxwell.R
\name{TI.MLE}
\alias{TI.MLE}
\title{Two-sided or equal-tailed tolerance Interval}
\usage{
TI.MLE(x, nr, p, cl, tails)
}
\arguments{
\item{x}{a numeric vercor of data input}

\item{nr}{numeric vector length one, number of simulation runs, defaults to 100,000.}

\item{p}{a numeric vector of probabilities}

\item{cl}{confidence level as in decimal, numeric value between 0 and 1.}

\item{tails}{logical; If TRUE, function calculate Equal-tailed tolerance factors. Otherwise, function calculate two-sided tolerance factor.}
}
\description{
This function computes (p, cl) two-sided or equal-tailed tolerance Interval
}
\examples{
#' # x<-c(135, 98, 114, 137, 138, 144, 99, 93, 115, 106, 132, 122, 94, 98, 127,
# 122, 102, 133, 114, 120, 93, 126, 119, 104, 119, 114, 125, 107, 98, 117, 111,
# 106, 108, 127, 126, 135, 112, 94, 127, 99, 120, 120, 121, 122, 96, 109, 123, 105)

For equal-tailed tolerance interval:
# TI.MLE(x,10^5,.9,.95,T)
For two-sided tolerance interval:
# TI.MLE(x,10^5,.9,.95,F)

}
\keyword{TI.MLE}
