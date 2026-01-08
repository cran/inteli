\name{lipois}
\alias{lipois}
\title{Likelihood Interval for Poisson Mean of Single Sample from Poisson Distribution}
\description{
  Likelihood interval for Poisson mean of single sample assuming Poisson distribution. Likelihood function plot is also shown. For faster calculation speed, plot option '0' is recommended.
}
\arguments{
  \item{event}{number of events; a numeric vector}
  \item{unit.time}{unit time for Poisson mean; default is 1}
  \item{data.size}{sample size. If not specified, size is assumed to be 1. For LI being more precise, sample size must be considered.}
  \item{conf.level}{confidence level for the interval estimation}
  \item{eps}{machine epsilon value; default is 1e-08}
  \item{k}{cutoff value for likelihood interval. If not specified, F-test is used.}
}
\value{
  \item{demo}{basic demographics of data}
  \item{LI}{estimated interval of Poisson mean by LI}
  \item{CI}{estimated interval of Poisson mean by CI}
  \item{plot}{likelihood-related plots}
}
\author{Kyun-Seop Bae <k@acr.kr>, Minkyu Kim <mkim@acr.kr>}
\examples{
  lipois(8)
  lipois(event = 8, unit.time = 3)
  lipois(event = 8, unit.time = 3, data.size = 6)
}
