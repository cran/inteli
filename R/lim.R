#' @importFrom stats binom.test nlminb poisson.test qchisq qf qnorm qt t.test uniroot var var.test
#' @importFrom graphics abline legend par points
#' @export

lim = function (data, conf.level = 0.95, df = 1, k)
{
  x <- data[!is.na(data)]
  n0 <- length(x)

  if (!is.numeric(x) | sum(is.infinite(x) > 0) | sum(is.nan(x)) > 0 | n0 < 3 |
      length(unique(x)) == 1) stop("Check the input!")

  m0 <- mean(x)

  if (!missing(k)) {
    logk <- log(k)
  } else {
    logk <- n0/2*log(1 + qf(conf.level, 1, n0 - df)/(n0 - df))
    logk <- min(logk, log(2/(1 - conf.level)))
  }

  LL <- function(th) -n0/2 * (log(2 * pi * sum((x - th)^2)/n0) + 1)
  maxLL <- sapply(m0, LL)

  O1 <- function(th) maxLL - sapply(th, LL) - logk
  O2 <- function(th) sapply(th, LL)
  O3 <- function(th) maxLL - sapply(th, LL)

  mLL <- uniroot(O1, c(-1e+08, m0), tol=1e-10)$root
  mUL <- uniroot(O1, c(m0, 1e+06), tol=1e-10)$root

  mLLc <- m0 - qt(0.5 + conf.level/2, n0 - 1) * sqrt(var(x)/n0)
  mULc <- m0 + qt(0.5 + conf.level/2, n0 - 1) * sqrt(var(x)/n0)

  demo <- c("Sample Size" = n0, "Cutoff Value k" = exp(logk), "maxLL" = maxLL)
  LI <- c("PE" = m0, "lower" = mLL, "upper" = mUL, "width" = mUL - mLL)
  CI <- c("PE" = m0, "lower" = mLLc, "upper" = mULc, "width" = mULc - mLLc)

  dm <- seq(mLL - 2, mUL + 2, length.out = 1e3)
  O1plot <- function(th) {
    plot(th, O1(th), type = "l",
         xlab = "Mean Value",
         ylab = "maxLL - LL - logk",
         main = "Objective Function (O1 type)")
    abline(h = 0, col = "red")
    abline(v = m0, lty=2)
    legend("topright",
           legend = c(paste("PE = ", format(m0, digits = 2)), "Zero Line"),
           lty = c(2, 1),
           col = c("black", "red"))
  }
  O2plot <- function(th) {
    plot(th, O2(th), type = "l",
         xlab = "Mean Value",
         ylab = "LL",
         main = "Log Likelihood Function (O2 type)")
    abline(h = maxLL, col = "blue")
    abline(v = m0, lty=2)
    abline(h = maxLL - logk, col = "red")
    legend("bottomright",
           legend = c(paste("PE = ", format(m0, digits = 2)),
                      paste("maxLL = ", format(maxLL, digits = 4)),
                      paste("maxLL-logk = ", format(maxLL - logk, digits = 4))),
           lty = c(2, 1, 1),
           col = c("black", "blue", "red"))
  }
  O3plot <- function(th) {
    plot(th, O3(th), type = "l",
         xlab = "Mean Value",
         ylab = "maxLL - LL",
         main = "Log LRT (O3 type)")
    abline(h = logk, col = "red")
    abline(v = m0, lty = 2)
    legend("topright",
           legend=c(paste("PE = ", format(m0, digits = 2)),
                    paste("logk = ", format(logk, digits = 4))),
           lty=c(2, 1),
           col=c("black", "red"))
  }
  par(mfrow = c(2,2))
  O1plot(dm); O2plot(dm); O3plot(dm)
  par(mfrow = c(1,1))
  return(list(demo = demo, LI = LI, CI = CI))
}
