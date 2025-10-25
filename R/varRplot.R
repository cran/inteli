#' varRplot: Plot of Variance Ratio Estimate by Likelihood Method
#'
#' This function plots a graph of interval estimation for a two group variance ratio by LI method, either in the log-likelihood function or the normalized log-likelihood value.
#' @param num.data A numeric vector functioning as a sample data, in a numerator position.
#' @param denom.data A numeric vector functioning as a sample data, in a denominator position.
#' @param logLRT A function type to be plotted. A default value "FALSE" refers to the log-likelihood function plot, while "TRUE" refers to the normalized log-likelihood ratio plot, or maxLL-LL.
#' @param conf.level A confidence level for CI method.
#' @param df A degree of freedom for LI method in terms of the denominator degree of freedom of F-test, as (n-df) of LRT, where n is the sum of sample sizes of input datum. For a variance ratio estimation, it is suggested to be 2.4.
#' @param low.scale A scaling factor for plotting the minimum value of x-axis, or a parameter value. The plot starts from "PE/low.scale". 5 is a default.
#' @param up.scale A scaling factor for plotting the maximum value of x-axis, or a parameter value. The plot starts from "PE*up.scale". 5 is a default.
#' @param k A cutoff value for LI method. Unless specified, F-test is used.
#' @return Plotted graph, either in the log-likelihood function or the normalized log-likelihood value
#' @examples
#' x <- rnorm(20, 0, 1)
#' y <- rnorm(40, 0, 1)
#' varRplot(x, y, FALSE)
#' varRplot(x, y, TRUE)
#'
#' @importFrom stats qf uniroot qchisq
#' @importFrom graphics abline legend
#' @export
varRplot = function (num.data, denom.data, logLRT = FALSE, conf.level = 0.95,
                     df = 2.4, low.scale = 5, up.scale = 5, k)

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
  Obj1 <- function(r) {
    th <- (n1v1 + r*n2v2)/r/n0
    ln2pith <- log(2*pi*th)
    -(n1*(log(r) + ln2pith) + n1v1/r/th + n2*ln2pith + n2v2/th)/2
  }
  Obj2 <- function(r) {
    th <- (n1v1 + r*n2v2)/r/n0
    ln2pith <- log(2*pi*th)
    maxLL + (n1*(log(r) + ln2pith) + n1v1/r/th + n2*ln2pith + n2v2/th)/2
  }
  dth <- seq(R0/low.scale, up.scale*R0, length.out=1e3)
  if (logLRT == FALSE){
    plot(dth, Obj1(dth), type="l", xlab="variance ratio estimation",
         ylab="log-likelihood function value")
    abline(h=maxLL, col="red")
    abline(v=R0, lty=2)
    abline(h=maxLL-logk, col="blue")
    legend("bottomright",
           legend=c(paste("Point Estimate = ", format(R0, digits=2)),
                    paste("maxLL = ", format(maxLL, digits=4)),
                    paste("maxLL-logk = ", format(maxLL-logk, digits=4))),
           lty=c(2, 1, 1),
           col=c("black", "red", "blue"))
  } else if (logLRT == TRUE) {
    plot(dth, Obj2(dth), type="l", xlab="variance ratio estimation",
         ylab="delta (maxLL - LL) value")
    abline(h=logk, col="red")
    abline(v=R0, lty=2)
    legend("bottomright",
           legend=c(paste("Point Estimate = ", format(R0, digits=2)),
                    paste("logk = ", format(logk, digits=4))),
           lty=c(2, 1),
           col=c("black", "red"))
  }
}
