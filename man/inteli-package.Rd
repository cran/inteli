\name{inteli-package}
\alias{inteli-package}
\alias{inteli}
\docType{package}
\title{Interval Estimation by Likelihoodist (LI) Compared to Frequentist (CI)}
\description{Parameter estimation via likelihood interval (LI) compared to traditional method (CI). This is the expanded version for 'LBI'- and 'wnl'-package, formulated by Kyun-Seop Bae <k@acr.kr>.}
\details{
Currently used CI method has its limitation when the test statistics are asymmetrical (chi-square test, F-test) or the model functions are non-linear. It can be overcome by using the likelihood functions for the interval estimation. 'inteli' package now supports interval estimation for the variance and variance ratio, and their likelihood function plots. Further estimating methods (especially for non-linear model functions) are upcoming. 
}
\author{Minkyu Kim <mkim@acr.kr>}
\references{
1. Wilks SS. The Large-sample Distribution of the Likelihood Ratio for Testing Composite Hypotheses. Ann Math Stat. 1938;9(1):60-62.

2. Ruppert D, Cressie N, Carroll RJ. A Transformation/Weighting Model for Estimating Michaelis-Menten Parameters. Cornell University Technical Report 796. 1988.

3. Pawitan Y. In All Likelihood: Statistical Modelling and Inference Using Likelihood. 2001.

4. Lee MH, Bae KS. Likelihood interval for nonlinear regression. Transl Clin Pharmacol. 2023.
}
\keyword{Packages}
\examples{
x <- rnorm(20, 0, 1)
y <- rnorm(40, 0, 1)

varE(x)
varEplot(x, FALSE)
varEplot(x, TRUE)

varR(x, y)
varRplot(x, y, FALSE)
varRplot(x, y, TRUE)
}
