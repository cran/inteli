\name{qliv}
\alias{qliv}
\title{Testing Function for Variance Ratio of Two Samples from Normal Distribution}
\description{
  Testing function for variance ratio of two samples assuming normal distribution by LI method. O2 type likelihood function plot is also shown. Smaller sample is allocated to be in the nominator for more accurate and precise estimation.
}
\arguments{
  \item{set.a}{a numeric vector of data values; mandatory, pre-allocation}
  \item{set.b}{a numeric vector of data values; mandatory, pre-allocation}
  \item{conf.level}{confidence level for the interval estimation}
  \item{df}{degree of freedom for likelihood interval; default is 2.4}
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
  x <- rnorm(30, 0, 1)
  y <- rnorm(20, 0, 1)
  qliv(x, y) # x placed in the denominator
}
