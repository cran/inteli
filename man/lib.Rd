\name{lib}
\alias{lib}
\title{Likelihood Interval for Probability of Single Sample from Binomial Distribution}
\description{
  Likelihood interval for probability of single sample assuming binomial distribution. Likelihood function plot is also shown. For faster calculation speed, plot option '0' is recommended.
}
\arguments{
  \item{event}{number of events from observed data; a numeric vector}
  \item{total}{total number of trials from observed data; a numeric vector}
  \item{conf.level}{confidence level for the interval estimation}
  \item{eps}{machine epsilon value; default is 1e-08}
  \item{k}{cutoff value for likelihood interval. If not specified, F-test is used.}
}
\value{
  \item{demo}{basic demographics of data}
  \item{LI}{estimated interval of probability (y/n) by LI}
  \item{CI}{estimated interval of probability (y/n) by CI}
  \item{plot}{likelihood-related plots}
}
\author{Kyun-Seop Bae <k@acr.kr>, Minkyu Kim <mkim@acr.kr>}
\examples{
  lib(8, 10)
}
