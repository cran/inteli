#' @importFrom stats binom.test nlminb poisson.test qchisq qf qnorm qt t.test uniroot var var.test
#' @importFrom graphics abline legend par points
#' @export

lipois = function (event, unit.time = 1, data.size = 1,
                   plot = "all", conf.level = 0.95, eps = 1e-08, k)
{
  ifelse (length(event) == 1, Lam0 <- event/unit.time,
          "! Event count should be a natural number")

  ifelse (data.size == 1, n <- 1, n <- data.size)

  maxLL <- n * Lam0 * (log(Lam0) - 1)

  if (!missing(k)) {
    logk <- log(k)
  } else if (n == 1) {
    logk <- log(2/(1 - conf.level))
  } else {
    logk <- n/2 * log(1 + qf(conf.level, 1, n - 1)/(n - 1))
    logk <- min(logk, log(2/(1 - conf.level)))
  }

  O1 <- function(lam, cutoff) {
    LL <- ifelse(lam > 0, n * Lam0 * log(lam) - n * lam, 0)
    maxLL - LL - cutoff
  }
  O2 <- function(lam, cutoff) ifelse(lam > 0, n * Lam0 * log(lam) - n * lam, 0)
  O3 <- function(lam, cutoff) {
    LL <- ifelse(lam > 0, n * Lam0 * log(lam) - n * lam, 0)
    maxLL - LL
  }

  if (Lam0 > 0) {
    pLL <- uniroot(O1, c(eps, Lam0), cutoff = logk)$root
  } else {
    pLL <- 0
  }

  if (!is.finite(maxLL)) {
    pUL <- Inf
  } else {
    pUL <- uniroot(O1, c(Lam0, 1e+09), cutoff = logk)$root
  }

  CIp <- poisson.test(x = event, T = unit.time, conf.level = conf.level)
  demo <- c("Poisson Mean" = Lam0,
            "Data Size" = ifelse(data.size == 1, 1, data.size),
            "Cutoff Value k" = exp(logk), "maxLL" = maxLL)
  LI <- c("Point Estimate" = Lam0, "lower" = pLL, "upper" = pUL,
          "width" = pUL - pLL)
  CI <- c("Point Estimate" = unname(CIp$estimate),
          "lower" = CIp$conf.int[1], "upper" = CIp$conf.int[2],
          "width" = CIp$conf.int[2] - CIp$conf.int[1])
  plot <- c("Current plot setting is" = plot)

  dp <- seq(0, pUL + 2, length.out = 1e3)
  if (plot == "all") {
    O1plot <- function(p) {
      plot(p, O1(p, logk), type = "l",
           xlab = "Poisson Mean Value",
           ylab = "maxLL - LL - logk",
           main = "Objective Function (O1 type)")
      abline(h = 0, col = "red")
      abline(v = Lam0, lty=2)
      legend("topright",
             legend = c(paste("PE = ", format(Lam0, digits = 2)), "Zero Line"),
             lty = c(2, 1),
             col = c("black", "red"))
    }
    O2plot <- function(p) {
      plot(p, O2(p, logk), type = "l",
           xlab = "Poisson Mean Value",
           ylab = "LL",
           main = "Log Likelihood Function (O2 type)")
      abline(h = maxLL, col = "blue")
      abline(v = Lam0, lty=2)
      abline(h = maxLL - logk, col = "red")
      legend("bottomright",
             legend=c(paste("PE = ", format(Lam0, digits = 2)),
                      paste("maxLL = ", format(maxLL, digits = 4)),
                      paste("maxLL-logk = ", format(maxLL-logk, digits = 4))),
             lty=c(2, 1, 1),
             col=c("black", "blue", "red"))
    }
    O3plot <- function(p) {
      plot(p, O3(p, logk), type = "l",
           xlab = "Poisson Mean Value",
           ylab = "maxLL - LL",
           main = "Log LRT (O3 type)")
      abline(h = logk, col = "red")
      abline(v = Lam0, lty = 2)
      legend("topright",
             legend = c(paste("PE = ", format(Lam0, digits = 2)),
                      paste("logk = ", format(logk, digits = 4))),
             lty = c(2, 1),
             col = c("black", "red"))
    }
    par(mfrow = c(2,2))
    O1plot(dp); O2plot(dp); O3plot(dp)
    par(mfrow = c(1,1))
  } else if (plot == "OBJ" | plot == "O1" | plot == 1) {
    O1plot(dp)
  } else if (plot == "OFV" | plot == "O2" | plot == 2) {
    O2plot(dp)
  } else if (plot == "LRT" | plot == "O3" | plot == 3) {
    O3plot(dp)
  } else {}

  return(list(demo = demo, LI = LI, CI = CI, plot = plot))
}
