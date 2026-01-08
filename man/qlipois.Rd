\name{qlipois}
\alias{qlipois}
\title{Testing Function for Poisson Mean from Poisson Distribution}
\description{
  Testing function for Poisson mean of single sample assuming Poisson distribution by LI method. O2 type likelihood function plot is also shown.
}
\arguments{
  \item{data}{a numeric vector of data values; mandatory}
  \item{test.val}{a numeric vector of testing value for Poisson mean; mandatory}
  \item{conf.level}{confidence level for the interval estimation}
  \item{eps}{machine epsilon value; default is 1e-08}
  \item{k}{cutoff value for likelihood interval. If not specified, F-test is used.}
}
\value{
  \item{demo}{basic demographics of data}
  \item{interval}{estimated interval by LI and CI.}
  \item{stat}{statistical values, such as p-value and likelihood ratio}
  \item{verdict}{estimation result as a sentence}
  \item{plot}{log likelihood function plot (O2 type)}
}
\author{Kyun-Seop Bae <k@acr.kr>, Minkyu Kim <mkim@acr.kr>}
\examples{
  qlipois(rnorm(20, 1, 1), 1)
}
