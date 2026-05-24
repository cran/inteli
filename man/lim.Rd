\name{lim}
\alias{lim}
\title{Likelihood Interval for Mean of Single Sample from Normal Distribution}
\description{
  Likelihood interval for mean of single sample assuming normal distribution. Likelihood function plot is also shown. Likelihood function plot is also shown. For faster calculation speed, plot option '0' is recommended.
}
\arguments{
  \item{data}{a numeric vector of data values; mandatory}
  \item{conf.level}{confidence level for the interval estimation}
  \item{df}{degree of freedom for likelihood interval; default is 1}
  \item{k}{cutoff value for likelihood interval. If not specified, F-test is used.}
}
\value{
  \item{demo}{basic demographics of data}
  \item{LI}{estimated interval of mean by LI}
  \item{CI}{estimated interval of mean by CI}
  \item{plot}{likelihood-related plots}
}
\author{Kyun-Seop Bae <k@acr.kr>, Minkyu Kim <mkim@acr.kr>}
\examples{
  set.seed(2026)
  lim(rnorm(30, 0, 1))
}
