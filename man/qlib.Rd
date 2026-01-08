\name{qlib}
\alias{qlib}
\title{Testing Function for Probability from Binomial Distribution}
\description{
  Testing function for variance ratio of two samples assuming normal distribution by LI method. O2 type likelihood function plot is also shown.
}
\arguments{
  \item{event}{number of events from observed data; a numeric vector}
  \item{total}{total number of trials from observed data; a numeric vector}
  \item{test.val}{a numeric vector of testing value for probability}
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
  qlib(8, 13, 0.7)
  qlib(3, 10, 0.8)
}
