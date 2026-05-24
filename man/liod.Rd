\name{liod}
\alias{liod}
\title{Likelihood Interval for Odds Ratio of 2-by-2 Contingency Table.}
\description{
  Likelihood interval for odds ratio of 2-by-2 contingency table. Likelihood function plot is also shown. For faster calculation speed, plot option '0' is recommended.
}
\arguments{
  \item{exposure.O_event.O}{a numeric vector of exposure (+) and event (+), 'a' of 2-by-2 table; mandatory}
  \item{exposure.O_TOTAL}{a numeric vector of exposure (+), 'a' + 'b' of 2-by-2 table; mandatory}
  \item{exposure.X_event.O}{a numeric vector of exposure (-) and event (-), 'c' of 2-by-2 table; mandatory}
  \item{exposure.X_TOTAL}{a numeric vector of exposure (-), 'c' + 'd' of 2-by-2 table; mandatory}
  \item{conf.level}{confidence level for the interval estimation}
  \item{eps}{machine epsilon value; default is 1e-08}
  \item{k}{cutoff value for likelihood interval. If not specified, F-test is used.}
}
\value{
  \item{demo}{basic demographics of data}
  \item{contin_2x2}{2-by-2 contingency table}
  \item{LI}{estimated interval of odds ratio by LI}
  \item{CI}{estimated interval of odds ratio by CI}
  \item{plot}{likelihood-related plots}
}
\author{Kyun-Seop Bae <k@acr.kr>, Minkyu Kim <mkim@acr.kr>}
\examples{
  liod(21, 43, 16, 67)
}
