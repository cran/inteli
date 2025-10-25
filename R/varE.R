#' varE: Calculate Variance Estimate
#'
#' This function computes the interval estimation for a single group variance by both LI and CI method.
#' @param data A numeric vector functioning as a sample data.
#' @param conf.level A confidence level for CI method.
#' @param df A degree of freedom for LI method in terms of the denominator degree of freedom of F-test, as (n-df) of LRT, where n is the sample size of input data. For a variance estimation, it is suggested to be 1.2.
#' @param lower A lower bound of "uniroot" for lower limit (LL) calculation. 1e-08 is a default.
#' @param upper A upper bound of "uniroot" for upper limit (UL) calculation. 1e+06 is a default.
#' @param k A cutoff value for LI method. Unless specified, F-test is used.
#' @return Point Estimate (PE), lower limit/bound (LL/LB), upper limit/bound (UL/UB), width, sample size, cutoff value k and maximum log-likelihood function value are calculated.
#' @examples
#' x <- rnorm(20, 0, 1)
#' varE(x)
#'
#' y <- rnorm(40, 0, 1)
#' varE(y)
#'
#' @importFrom stats qf uniroot qchisq
#' @export
varE = function (data, conf.level = 0.95, df = 1.2,
                 lower = 1e-08, upper = 1e+06, k)
{
  x <- data[!is.na(data)]
  n0 <- length(x)
  if (!is.numeric(x) | sum(is.infinite(x) > 0) | sum(is.nan(x)) >
      0 | n0 < 3 | length(unique(x)) == 1)
    stop("Check the input!")
  n0v0 <- sum((x - mean(x))^2)
  v0 <- n0v0/n0
  v0c <- n0v0/(n0-1)
  maxLL <- -n0 * (log(2 * pi * v0) + 1)/2
  if (!missing(k)) {logk <- log(k)}
  else {logk <- n0/2 * log(1 + qf(conf.level, 1, n0 - df)/(n0 - df))
  logk <- min(logk, log(2/(1 - conf.level)))}
  O2 <- function(th) maxLL + (n0 * log(2 * pi * th) + n0v0/th)/2 - logk
  varLL <- uniroot(O2, c(lower, v0))$root
  varUL <- uniroot(O2, c(v0, upper * v0))$root
  varLB <- n0v0/qchisq(0.5 + conf.level/2, n0 - 1)
  varUB <- n0v0/qchisq(0.5 - conf.level/2, n0 - 1)
  LI <- c("Point Estimate" = v0, "lower" = varLL, "upper" = varUL,
          "width" = varUL - varLL)
  CI <- c("Point Estimate" = v0c, "lower" = varLB, "upper" = varUB,
          "width" = varUB - varLB)
  size <- c("Number of Data" = n0, "df for CI" = n0-1)
  cutoff <- c("Cutoff Value k" = exp(logk), "Cutoff Value logk" = logk)
  maxLL <- c("Maximum log-Likelihood Fx Value" = maxLL)
  z <- list(LI = LI, CI = CI, size = size, cutoff = cutoff, maxLL = maxLL)
  z
}
