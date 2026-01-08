#' @importFrom stats binom.test nlminb poisson.test qchisq qf qnorm qt t.test uniroot var var.test
#' @importFrom graphics abline legend par points
#' @export

liv = function (data, plot = "all", conf.level = 0.95, df = 1.2, k)
{
  x <- data[!is.na(data)]
  n0 <- length(x)

  if (!is.numeric(x) | sum(is.infinite(x) > 0) | sum(is.nan(x)) >
      0 | n0 < 3 | length(unique(x)) == 1)
    stop("Check the input!")

  n0v0 <- sum((x - mean(x))^2)
  v0 <- n0v0 / n0
  v0c <- n0v0 / (n0 - 1)

  if (!missing(k)) {
    logk <- log(k)
    } else {
      logk <- n0/2 * log(1 + qf(conf.level, 1, n0 - df)/(n0 - df))
      logk <- min(logk, log(2/(1 - conf.level)))
    }

  O2 <- function(th) -(n0 * log(2 * pi * th) + n0v0/th)/2
  maxLL <- O2(v0)
  O1 <- function(th) maxLL - O2(th) - logk
  O3 <- function(th) maxLL - O2(th)

  varLL <- uniroot(O1, c(1e-08, v0))$root
  varUL <- uniroot(O1, c(v0, 1e+06))$root

  varLB <- n0v0 / qchisq(0.5 + conf.level/2, n0 - 1)
  varUB <- n0v0 / qchisq(0.5 - conf.level/2, n0 - 1)

  demo <- c("Sample Size" = n0, "Cutoff Value k" = exp(logk), "maxLL" = maxLL)
  LI <- c("Point Estimate" = v0, "lower" = varLL, "upper" = varUL,
          "width" = varUL - varLL)
  LI.sd <- c("Point Estimate" = sqrt(v0), "lower" = sqrt(varLL),
             "upper" = sqrt(varUL))
  CI <- c("Point Estimate" = v0c, "lower" = varLB, "upper" = varUB,
          "width" = varUB - varLB)
  plot <- c("Current plot setting is" = plot)

  dth <- seq(varLL/2, varUL * 2, length.out = 1e3)
  if (plot == "all") {
    O1plot <- function(th) {
      plot(th, O1(th), type = "l",
           xlab = "Variance Value",
           ylab = "maxLL - LL - logk",
           main = "Objective Function (O1 type)")
      abline(h = 0, col = "red")
      abline(v = v0, lty=2)
      legend("topright",
             legend = c(paste("PE = ", format(v0, digits = 2)), "Zero Line"),
             lty = c(2, 1),
             col = c("black", "red"))
    }
    O2plot <- function(th) {
      plot(th, O2(th), type = "l",
           xlab = "Variance Value",
           ylab = "LL",
           main = "Log Likelihood Function (O2 type)")
      abline(h = maxLL, col = "blue")
      abline(v = v0, lty=2)
      abline(h = maxLL - logk, col = "red")
      legend("bottomright",
             legend = c(paste("PE = ", format(v0, digits = 2)),
                        paste("maxLL = ", format(maxLL, digits = 4)),
                        paste("maxLL-logk = ", format(maxLL - logk, digits = 4))),
             lty = c(2, 1, 1),
             col = c("black", "blue", "red"))
    }
    O3plot <- function(th) {
      plot(th, O3(th), type = "l",
           xlab = "Variance Value",
           ylab = "maxLL - LL",
           main = "Log LRT (O3 type)")
      abline(h = logk, col = "red")
      abline(v = v0, lty = 2)
      legend("topright",
             legend=c(paste("PE = ", format(v0, digits = 2)),
                      paste("logk = ", format(logk, digits = 4))),
             lty=c(2, 1),
             col=c("black", "red"))
    }
    par(mfrow = c(2,2))
    O1plot(dth); O2plot(dth); O3plot(dth)
    par(mfrow = c(1,1))
  } else if (plot == "OBJ" | plot == "O1" | plot == 1) {
    O1plot(dth)
  } else if (plot == "OFV" | plot == "O2" | plot == 2) {
    O2plot(dth)
  } else if (plot == "LRT" | plot == "O3" | plot == 3) {
    O3plot(dth)
  } else {}

  return(list(demo = demo, LI = LI, LI.sd = LI.sd, CI = CI, plot = plot))
}
