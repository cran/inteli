\name{inteli-package}
\alias{inteli-package}
\alias{inteli}
\docType{package}
\title{Interval Estimation by Likelihood Method}
\description{
Parameter estimation via likelihood interval (LI) compared to conventional method (CI).
}
\details{
  Currently used CI method has its limitation when the test statistics are asymmetrical (chi-square test, F-test) or the model functions are non-linear. It can be overcome by using the likelihood functions for the interval estimation. 'inteli' package now supports interval estimation for the mean, variance, variance ratio, binomial distribution, Poisson distribution, odds ratio, risk difference, relative risk and their likelihood function plots. Testing functions are also provided.
}
\author{
\strong{Maintainer}: Minkyu Kim \email{mkim@acr.kr}

Authors:
\itemize{
  \item Kyun-Seop Bae \email{k@acr.kr} (\href{https://orcid.org/0000-0001-7399-5879}{ORCID})
  \item Minkyu Kim \email{mkim@acr.kr}
}
}
\references{
\enumerate{
\item Wilks SS. The Large-sample Distribution of the Likelihood Ratio for Testing Composite Hypotheses. Ann Math Stat. 1938;9(1):60-62.
\item Edwards AWF. Likelihood. 1972.
\item Fisher RA. Statistical Methods and Scientific Inference. 3e. 1973.
\item Bates DM, Watts DG. Nonlinear Regression Analysis and its Application. 1988.
\item Ruppert D, Cressie N, Carroll RJ. A Transformation/Weighting Model for Estimating Michaelis-Menten Parameters. Cornell University Technical Report 796. 1988.
\item Royall R. Statistical Evidence. 1997.
\item Pinheiro JC, Bates DM. Mixed Effects Models in S and S-PLUS. 2000.
\item Pawitan Y. In All Likelihood: Statistical Modelling and Inference Using Likelihood. 2001.
\item Lehmann EL. Fisher, Nayman, and the Creation of Classical Statistics. 2011.
\item Rohde CA. Introductory Statistical Inference with the Likelihood Function. 2014.
\item Held L, Bove DS. Likelihood and Bayesian Inference. 2020.
\item Lee MH, Bae KS. Likelihood interval for nonlinear regression. 2023.
}
}
