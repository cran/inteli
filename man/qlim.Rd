\name{qlim}
\alias{qlim}
\title{Testing Function for Mean Difference of Two Samples from Normal Distribution by LI Method}
\description{
  Testing function for mean difference of two samples assuming normal distribution by LI method. O2 type likelihood function plot is also shown. 
}
\arguments{
  \item{set.a}{a numeric vector of former data values; mandatory}
  \item{set.b}{a numeric vector of latter data values; mandatory}
  \item{conf.level}{confidence level for the interval estimation}
  \item{df}{degree of freedom for likelihood interval; default is 1}
  \item{k}{cutoff value for likelihood interval. If not specified, F-test is used.}
}
\value{
  \item{demo}{basic demographics of data}
  \item{interval}{estimated interval by LI and CI. PE is mean(set.a) - mean(set.b).}
  \item{stat}{statistical values, such as p-value and likelihood ratio}
  \item{verdict}{estimation result as a sentence}
  \item{plot}{log likelihood function plot (O2 type)}
}
\author{Kyun-Seop Bae <k@acr.kr>, Minkyu Kim <mkim@acr.kr>}
\examples{
  x <- rnorm(30, 0, 1)
  y <- rnorm(20, 0, 1)
  qlim(x, y)
}
