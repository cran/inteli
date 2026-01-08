\name{lird}
\alias{lird}
\title{Likelihood Interval for Risk Difference of 2-by-2 Contingency Table.}
\description{
  Likelihood interval for risk difference of 2-by-2 contingency table. Likelihood function plot is also shown. For faster calculation speed, plot option '0' is recommended.
}
\arguments{
  \item{exposure.O_event.O}{a numeric vector of exposure (+) and event (+), 'a' of 2-by-2 table; mandatory}
  \item{exposure.O_TOTAL}{a numeric vector of exposure (+), 'a' + 'b' of 2-by-2 table; mandatory}
  \item{exposure.X_event.O}{a numeric vector of exposure (-) and event (-), 'c' of 2-by-2 table; mandatory}
  \item{exposure.X_TOTAL}{a numeric vector of exposure (-), 'c' + 'd' of 2-by-2 table; mandatory}
  \item{plot}{likelihood function plots by LI; default is "all". 'OBJ', 'O1', '1' is objective function plot; 'OFV', 'O2', '2' is log likelihood function plot; 'LRT', 'O3', '3' is log LRT function plot.}
  \item{conf.level}{confidence level for the interval estimation}
  \item{eps}{machine epsilon value; default is 1e-08}
  \item{k}{cutoff value for likelihood interval. If not specified, F-test is used.}
}
\value{
  \item{demo}{basic demographics of data}
  \item{contin_2x2}{2-by-2 contingency table}
  \item{LI}{estimated interval of risk difference by LI}
  \item{CI}{estimated interval of risk difference by CI}
  \item{plot}{likelihood-related plots}
}
\author{Kyun-Seop Bae <k@acr.kr>, Minkyu Kim <mkim@acr.kr>}
\examples{
  lird(21, 43, 16, 67)
  lird(21, 43, 16, 67, plot = 0) # faster calculation
}
