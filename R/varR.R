#' varR: Calculate Variance Ratio Estimate
#'
#' This function computes the interval estimation for a two group variance ratio by both LI and CI method.
#' @param num.data A numeric vector functioning as a sample data, in a numerator position.
#' @param denom.data A numeric vector functioning as a sample data, in a denominator position.
#' @param conf.level A confidence level for CI method.
#' @param df A degree of freedom for LI method in terms of the denominator degree of freedom of F-test, as (n-df) of LRT, where n is the sum of sample sizes of input datum. For a variance ratio estimation, it is suggested to be 2.4.
#' @param lower A lower bound of "uniroot" for lower limit (LL) calculation. 1e-08 is a default.
#' @param upper A upper bound of "uniroot" for upper limit (UL) calculation. 1e+06 is a default.
#' @param k A cutoff value for LI method. Unless specified, F-test is used.
#' @return Point Estimate (PE), lower limit/bound (LL/LB), upper limit/bound (UL/UB), width, sample size, cutoff value k and maximum log-likelihood function value are calculated.
#' @examples
#' x <- rnorm(20, 0, 1)
#' y <- rnorm(40, 0, 1)
#' varR(x, y)
#'
#' @importFrom stats qf uniroot qchisq
#' @export
varR = function (num.data, denom.data, conf.level = 0.95, df = 2.4,
                 lower = 1e-08, upper = 1e+06, k)
{
  x <- num.data[!is.na(num.data)]
  y <- denom.data[!is.na(denom.data)]
  n1 <- length(x)
  n2 <- length(y)
  n1v1 <- sum((x - mean(x))^2)
  n2v2 <- sum((y - mean(y))^2)
  v1 <- n1v1/n1
  v2 <- n2v2/n2
  v1c <- n1v1/(n1-1)
  v2c <- n2v2/(n2-1)
  if (!is.numeric(x) | sum(is.infinite(x) > 0) | sum(is.nan(x)) >
      0 | n1 < 3 | length(unique(x)) == 1)
    stop("Check the input!")
  if (!is.numeric(y) | sum(is.infinite(y) > 0) | sum(is.nan(y)) >
      0 | n2 < 3 | length(unique(y)) == 1)
    stop("Check the input!")
  R0 <- v1/v2
  R0c <- v1c/v2c
  n0 <- n1 + n2
  maxLL <- -(n1 * (log(2 * pi * v1) + 1) + n2 * (log(2 * pi * v2) + 1))/2
  if (!missing(k)) {logk <- log(k)}
  else {logk <- n0/2 * log(1 + qf(conf.level, 1, n0 - df)/(n0 - df))
  logk <- min(logk, log(2/(1 - conf.level)))}
  Obj <- function(r) {
    th <- (n1v1 + r*n2v2)/r/n0
    ln2pith <- log(2*pi*th)
    mLL <- -(n1*(log(r) + ln2pith) + n1v1/r/th + n2*ln2pith + n2v2/th)/2
    maxLL - mLL - logk
  }
  varLL <- uniroot(Obj, c(lower, R0))$root
  varUL <- uniroot(Obj, c(R0, upper))$root
  varLB <- R0c/qf(0.5 + conf.level/2, n1-1, n2-1)
  varUB <- R0c/qf(0.5 - conf.level/2, n1-1, n2-1)
  LI <- c("Point Estimate" = R0, "lower" = varLL, "upper" = varUL,
          "width" = varUL - varLL)
  CI <- c("Point Estimate" = R0c, "lower" = varLB, "upper" = varUB,
          "width" = varUB - varLB)
  size <- c("Number of Num. Data" = n1, "Number of Denom. Data" = n2,
            "Num. df for CI" = n1-1, "Denom. df for CI" = n2-1)
  cutoff <- c("Cutoff Value k" = exp(logk), "Cutoff Value logk" = logk)
  maxLL <- c("Maximum log-Likelihood Fx Value" = maxLL)
  z <- list(LI = LI, CI = CI, size = size, cutoff = cutoff, maxLL = maxLL)
  z
}
