\name{lir}
\alias{lir}
\title{Likelihood Interval for Variance Ratio of Two Samples from Normal Distribution}
\description{
  Likelihood interval for variance ratio of two samples assuming normal distribution. Likelihood function plot is also shown. For faster calculation speed, plot option '0' is recommended. Sample size is not allocated.
}
\arguments{
  \item{num.data}{a numeric vector of data values, variance value at the numerator position; mandatory}
  \item{denom.data}{a numeric vector of data values, variance value at the denominator position; mandatory}
  \item{plot}{likelihood function plots by LI; default is "all". 'OBJ', 'O1', '1' is objective function plot; 'OFV', 'O2', '2' is log likelihood function plot; 'LRT', 'O3', '3' is log LRT function plot.}
  \item{conf.level}{confidence level for the interval estimation}
  \item{df}{degree of freedom for likelihood interval; default is 2.4}
  \item{k}{cutoff value for likelihood interval. If not specified, F-test is used.}
}
\value{
  \item{demo}{basic demographics of data}
  \item{LI}{estimated interval of variance ratio by LI}
  \item{LI.sdR}{estimated interval of standard deviation ratio by LI}
  \item{CI}{estimated interval of variance by CI}
  \item{plot}{likelihood-related plots}
}
\author{Kyun-Seop Bae <k@acr.kr>, Minkyu Kim <mkim@acr.kr>}
\examples{
  lir(rnorm(30, 0, 1), rnorm(20, 0, 1))
}
