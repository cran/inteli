#' @importFrom stats binom.test nlminb poisson.test qchisq qf qnorm qt t.test uniroot var var.test
#' @importFrom graphics abline legend par points
#' @export

lib = function (event, total, plot = "all", conf.level = 0.95, eps = 1e-08, k)
{
  y <- event
  n <- total

  if (length(y) != 1 | length(n) != 1)
    stop("Event counts should be a sinlge natural number!")
  if (any(y < 0) | any(n < 0) | any(!is.finite(y)) | any(!is.finite(n)))
    stop("Check the event counts and(or) total counts!")
  if (!missing(k)) {
    logk <- log(k)
  } else if (n == 1) {
    logk <- log(2/(1 - conf.level))
  } else {
    logk <- n/2 * log(1 + qf(conf.level, 1, n - 1)/(n - 1))
    logk <- min(logk, log(2/(1 - conf.level)))
  }

  p0 <- y/n

  if (y == 0 | y == n) {
    LE <- logk/sqrt(n)
    maxLL <- 0
  } else {
    LE <- logk * sqrt(p0 * (1 - p0))
    maxLL <- y * log(p0) + (n - y) * log(1 - p0)
  }

  O1 <- function(p) {
    LL <- ifelse(p < 1e-300 | p > 1 - 1e-300, 0, y * log(p) +
                   (n - y) * log(1 - p))
    maxLL - LL - logk
  }
  O2 <- function(p) {
    ifelse(p < 1e-300 | p > 1 - 1e-300, 0, y * log(p) + (n - y) * log(1 - p))
  }
  O3 <- function(p) {
    LL <- ifelse(p < 1e-300 | p > 1 - 1e-300, 0, y * log(p) +
                  (n - y) * log(1 - p))
    maxLL - LL
  }

  if (p0 < eps) {
    bLL <- 0
  } else {
    rTemp <- try(uniroot(O1, c(max(eps, p0 - LE), p0 + eps)), silent = TRUE)
    if (!inherits(rTemp, "try-error")) {
      bLL <- rTemp$root
    } else {
      bLL <- max(p0 - LE, 0)
    }
  }

  if (p0 > 1 - eps) {
    bUL <- 1
  } else {
    rTemp = try(uniroot(O1, c(p0 - eps, min(p0 + LE, 1 - eps))), silent = TRUE)
    if (!inherits(rTemp, "try-error")) {
      bUL = rTemp$root
    } else {
      bUL <- min(p0 + LE, 1)
    }
  }

  CIb <- binom.test(y, n)

  demo <- c("Event Counts" = y, "Total Counts" = n,
            "Cutoff Value k" = exp(logk), "maxLL" = maxLL)
  LI <- c("Point Estimate" = p0, "lower" = bLL, "upper" = bUL,
          "width" = bUL - bLL)
  CI <- c("Point Estimate" = unname(CIb$estimate),
          "lower" = CIb$conf.int[1], "upper" = CIb$conf.int[2],
          "width" = CIb$conf.int[2] - CIb$conf.int[1])
  plot <- c("Current plot setting is" = plot)

  dp <- seq(0, 1, length.out = 1e3)
  if (plot == "all") {
    O1plot <- function(p) {
      plot(p, O1(p), type = "l",
           xlab = "Proportion Value",
           ylab = "maxLL - LL - logk",
           main = "Objective Function (O1 type)")
      abline(h = 0, col = "red")
      abline(v = p0, lty = 2)
      legend("topright",
             legend= c (paste("PE = ", format(p0, digits = 2)), "Zero Line"),
             lty = c(2, 1),
             col = c("black", "red"))
    }
    O2plot <- function(p) {
      plot(p, O2(p), type = "l",
           xlab = "Proportion Value",
           ylab = "LL",
           main = "Log Likelihood Function (O2 type)")
      abline(h = maxLL, col = "blue")
      abline(v = p0, lty = 2)
      abline(h = maxLL - logk, col = "red")
      legend("bottomright",
             legend = c(paste("PE = ", format(p0, digits = 2)),
                        paste("maxLL = ", format(maxLL, digits = 4)),
                        paste("maxLL-logk = ", format(maxLL-logk, digits = 4))),
             lty = c(2, 1, 1),
             col = c("black", "blue", "red"))
    }
    O3plot <- function(p) {
      plot(p, O3(p), type = "l",
           xlab = "Proportion Value",
           ylab = "maxLL - LL",
           main = "Log LRT (O3 type)")
      abline(h = logk, col = "red")
      abline(v = p0, lty = 2)
      legend("topright",
             legend = c(paste("PE = ", format(p0, digits = 2)),
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
