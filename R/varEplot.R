#' varEplot: Plot of Variance Estimate by Likelihood Method
#'
#' This function plots a graph of interval estimation for a single group variance by LI method, either in the log-likelihood function or the normalized log-likelihood value.
#' @param data A numeric vector functioning as a sample data.
#' @param logLRT A function type to be plotted. A default value "FALSE" refers to the log-likelihood function plot, while "TRUE" refers to the normalized log-likelihood ratio plot, or maxLL-LL.
#' @param conf.level A confidence level for CI method.
#' @param df A degree of freedom for LI method in terms of the denominator degree of freedom of F-test, as (n-df) of LRT, where n is the sample size of input data. For a variance estimation, it is suggested to be 1.2.
#' @param low.scale A scaling factor for plotting the minimum value of x-axis, or a parameter value. The plot starts from "PE/low.scale". 3 is a default.
#' @param up.scale A scaling factor for plotting the maximum value of x-axis, or a parameter value. The plot starts from "PE*up.scale". 5 is a default.
#' @param k A cutoff value for LI method. Unless specified, F-test is used.
#' @return Plotted graph, either in the log-likelihood function or the normalized log-likelihood value.
#' @examples
#' x <- rnorm(20, 0, 1)
#' varEplot(x, FALSE)
#' varEplot(x, TRUE)
#'
#' y <- rnorm(40, 0, 1)
#' varEplot(y, FALSE)
#' varEplot(y, TRUE)
#'
#' @importFrom stats qf uniroot qchisq
#' @importFrom graphics abline legend
#' @export
varEplot = function (data, logLRT = FALSE, conf.level = 0.95, df = 1.2,
                     low.scale = 3, up.scale = 5, k)
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
  O1 <- function(th) - (n0 * log(2 * pi * th) + n0v0/th)/2
  O2 <- function(th) maxLL + (n0 * log(2 * pi * th) + n0v0/th)/2
  dth <- seq(v0/low.scale, up.scale*v0, length.out=1e3)
  if (logLRT == FALSE){
  plot(dth, O1(dth), type="l", xlab="variance estimation",
       ylab="log-likelihood function value")
  abline(h=maxLL, col="red")
  abline(v=v0, lty=2)
  abline(h=maxLL-logk, col="blue")
  legend("bottomright",
         legend=c(paste("Point Estimate = ", format(v0, digits=2)),
                  paste("maxLL = ", format(maxLL, digits=4)),
                  paste("maxLL-logk = ", format(maxLL-logk, digits=4))),
         lty=c(2, 1, 1),
         col=c("black", "red", "blue"))
  } else if (logLRT == TRUE) {
    plot(dth, O2(dth), type="l", xlab="variance estimation",
         ylab="delta (maxLL - LL) value")
    abline(h=logk, col="red")
    abline(v=v0, lty=2)
    legend("bottomright",
           legend=c(paste("Point Estimate = ", format(v0, digits=2)),
                    paste("logk = ", format(logk, digits=4))),
           lty=c(2, 1),
           col=c("black", "red"))
  }
}
