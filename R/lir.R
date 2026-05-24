#' @importFrom stats binom.test nlminb poisson.test qchisq qf qnorm qt t.test uniroot var var.test
#' @importFrom graphics abline legend par points
#' @export

lir = function (num.data, denom.data, conf.level = 0.95, df.t = 2.4, k)
{
  x <- num.data[!is.na(num.data)]
  y <- denom.data[!is.na(denom.data)]
  n1 <- length(x)
  n2 <- length(y)

  if (!is.numeric(x) | sum(is.infinite(x) > 0) | sum(is.nan(x)) >
      0 | n1 < 3 | length(unique(x)) == 1)
    stop("Check the input!")
  if (!is.numeric(y) | sum(is.infinite(y) > 0) | sum(is.nan(y)) >
      0 | n2 < 3 | length(unique(y)) == 1)
    stop("Check the input!")

  n1v1 <- sum((x - mean(x))^2)
  n2v2 <- sum((y - mean(y))^2)
  v1 <- n1v1 / n1
  v2 <- n2v2 / n2
  v1c <- n1v1 / (n1 - 1)
  v2c <- n2v2 / (n2 - 1)
  R0 <- v1 / v2
  R0c <- v1c / v2c
  n0 <- n1 + n2

  if (n1 < 30 & n2 < 30) {
    df <- max(df.t, 1.75 * max(n1, n2)/min(n1, n2))
  } else df <- max(df.t, 1.5 * max(n1, n2)/min(n1, n2))

  if (!missing(k)) {
    logk <- log(k)
    } else {
      logk <- n0 / 2 * log(1 + qf(conf.level, 1, n0 - df)/(n0 - df))
      logk <- min(logk, log(2/(1 - conf.level)))
    }

  O2 <- function(r) {
    th <- (n1v1 + r * n2v2)/r/n0
    ln2pith <- log(2 * pi * th)
    -(n1 * (log(r) + ln2pith) + n1v1/r/th + n2 * ln2pith + n2v2/th)/2
  }
  maxLL <- -(n1 * (log(2 * pi * v1) + 1) + n2 * (log(2 * pi * v2) + 1))/2
  O1 <- function(r) maxLL - O2(r) - logk
  O3 <- function(r) maxLL - O2(r)

  varLL <- uniroot(O1, c(1e-08, R0))$root
  varUL <- uniroot(O1, c(R0, 1e+06))$root

  varLB <- R0c / qf(0.5 + conf.level/2, n1 - 1, n2 - 1)
  varUB <- R0c / qf(0.5 - conf.level/2, n1 - 1, n2 - 1)

  demo <- c("Num. Size" = n1, "Denom. Size" = n2,
            "Cutoff Value k" = exp(logk), "maxLL" = maxLL)
  LI <- c("Point Estimate" = R0, "lower" = varLL, "upper" = varUL,
          "width" = varUL - varLL)
  LI.sdR <- c("Point Estimate" = sqrt(R0), "lower" = sqrt(varLL),
              "upper" = sqrt(varUL))
  CI <- c("Point Estimate" = R0c, "lower" = varLB, "upper" = varUB,
          "width" = varUB - varLB)

  dr <- seq(varLL/2, varUL * 2, length.out = 1e3)

  O1plot <- function(r) {
    plot(r, O1(r), type = "l",
         xlab = "Variance Ratio Value",
         ylab = "maxLL - LL - logk",
         main = "Objective Function (O1 type)")
    abline(h = 0, col = "red")
    abline(v = R0, lty=2)
    legend("topright",
           legend = c(paste("PE = ", format(R0, digits = 2)), "Zero Line"),
           lty = c(2, 1),
           col = c("black", "red"))
  }
  O2plot <- function(r) {
    plot(r, O2(r), type = "l",
         xlab = "Variance Ratio Value",
         ylab = "LL",
         main = "Log Likelihood Function (O2 type)")
    abline(h = maxLL, col = "blue")
    abline(v = R0, lty=2)
    abline(h = maxLL - logk, col = "red")
    legend("bottomright",
           legend = c(paste("PE = ", format(R0, digits=2)),
                      paste("maxLL = ", format(maxLL, digits=4)),
                      paste("maxLL-logk = ", format(maxLL-logk, digits=4))),
           lty = c(2, 1, 1),
           col = c("black", "blue", "red"))
  }
  O3plot <- function(r) {
    plot(r, O3(r), type = "l",
         xlab = "Variance Ratio Value",
         ylab = "maxLL - LL",
         main = "Log LRT (O3 type)")
    abline(h = logk, col = "red")
    abline(v = R0, lty = 2)
    legend("topright",
           legend = c(paste("PE = ", format(R0, digits = 2)),
                      paste("logk = ", format(logk, digits = 4))),
           lty = c(2, 1),
           col = c("black", "red"))
  }
  par(mfrow = c(2,2))
  O1plot(dr); O2plot(dr); O3plot(dr)
  par(mfrow = c(1,1))
  return(list(demo = demo, LI = LI, LI.sdR = LI.sdR, CI = CI))
}
