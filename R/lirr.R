#' @importFrom stats binom.test nlminb poisson.test qchisq qf qnorm qt t.test uniroot var var.test
#' @importFrom graphics abline legend par points
#' @export

lirr = function (exposure.O_event.O, exposure.O_TOTAL,
                 exposure.X_event.O, exposure.X_TOTAL,
                 plot = "all", conf.level = 0.95, eps = 1e-08, k)
{
  y1 <- exposure.O_event.O
  n1 <- exposure.O_TOTAL
  y2 <- exposure.X_event.O
  n2 <- exposure.X_TOTAL

  if (length(y1) > 1)
    stop("Requires sinlge strata!")
  if (any(c(y1, n1 - y1, y2, n2 - y2) < 0) | n1 * n2 == 0)
    stop("Check the input!")

  p1 <- y1/n1
  p2 <- y2/n2
  d1 <- n1 - y1
  d2 <- n2 - y2
  RR0 <- p1/p2
  n <- n1 + n2

  if (!missing(k)) {
    logk <- log(k)
  } else if (n == 1) {
    logk <- log(2/(1 - conf.level))
  } else {
    logk <- n/2 * log(1 + qf(conf.level, 1, n - 1)/(n - 1))
    logk <- min(logk, log(2/(1 - conf.level)))
  }

  LE <- 1e4 * (exp(logk)/(n1 + n2))
  p1 <- min(1 - eps, max(eps, p1))
  p2 <- min(1 - eps, max(eps, p2))

  ifelse (y1 == 0 | y1 == n1 | y2 == 0 | y2 == n2, maxLL <- 0,
          maxLL <- y1 * log(p1) + d1 * log(1 - p1) + y2 * log(p2) + d2 * log(1 - p2))

  O1 <- function(rr) {
    A <- n * rr
    B <- n1 * rr + y1 + n2 + y2 * rr
    C1 <- y1 + y2
    p2t <- min(1 - eps, max(eps, (B - sqrt(B * B - 4 * A * C1))/(2 * A)))
    p1t <- min(1 - eps, max(eps, p2t * rr))
    LL <- y1 * log(p1t) + d1 * log(1 - p1t) + y2 * log(p2t) + d2 * log(1 - p2t)
    maxLL - LL - logk
  }

  LLt <- function(rr) {
    A <- n * rr
    B <- n1 * rr + y1 + n2 + y2 * rr
    C1 <- y1 + y2
    p2t <- (B - sqrt(B * B - 4 * A * C1))/(2 * A)
    p1t <- p2t * rr
    y1 * log(p1t) + d1 * log(1 - p1t) + y2 * log(p2t) + d2 * log(1 - p2t)
  }

  O1t <- function(rr) maxLL - LLt(rr) - logk
  O2t <- function(rr) LLt(rr)
  O3t <- function(rr) maxLL - LLt(rr)

  if (y1 == 0 & y2 == 0) {
    rLL <- 0
    rUL <- Inf
  } else if (y1 == 0 & y2 > 0 & y2 < n2) {
    rLL <- 0
    rTemp <- try(uniroot(O1, c(eps, 1e9)), silent = T)
    ifelse (!inherits(rTemp, "try-error"), rUL <- rTemp$root, rUL <- 1e+09)
  } else if (y1 == 0 & y2 == n2) {
    rLL <- 0
    rTemp <- try(uniroot(O1, c(eps, 1e9)), silent = T)
    ifelse (!inherits(rTemp, "try-error"), rUL <- rTemp$root, rUL <- Inf)
  } else if (y1 > 0 & y1 < n1 & y2 == 0) {
    rTemp <- try(uniroot(O1, c(eps, 100)), silent = T)
    ifelse (!inherits(rTemp, "try-error"), rLL <- rTemp$root, rLL <- 0)
    rUL <- Inf
  } else if (y1 > 0 & y1 < n1 & y2 > 0 & y2 < n2) {
    rTemp <- try(uniroot(O1, c(eps, RR0 + eps)), silent = T)
    ifelse (!inherits(rTemp, "try-error"),  rLL <- rTemp$root, rLL <- 0)
    rTemp <- try(uniroot(O1, c(RR0 - eps, RR0 + LE)), silent = T)
    ifelse (!inherits(rTemp, "try-error"), rUL <- rTemp$root, rUL <- Inf)
  } else if (y1 > 0 & y1 < n1 & y2 == n2) {
    rTemp <- try(uniroot(O1, c(eps, RR0 + eps)), silent = T)
    ifelse (!inherits(rTemp, "try-error"), rLL <- rTemp$root, rLL <- 0)
    rTemp <- try(uniroot(O1, c(RR0 - eps, 1e9)), silent = T)
    ifelse (!inherits(rTemp, "try-error"), rUL <- rTemp$root, rUL <- Inf)
  } else if (y1 == n1 & y2 == 0) {
    rTemp <- try(uniroot(O1, c(eps, 1e9)), silent = T)
    ifelse (!inherits(rTemp, "try-error"), rLL <- 1/rTemp$root, rLL <- 0)
    rUL <- Inf
  } else if (y1 == n1 & y2 > 0 & y2 < n2) {
    rTemp <- try(uniroot(O1, c(eps, RR0 + eps)), silent = T)
    ifelse (!inherits(rTemp, "try-error"), rLL <- rTemp$root, rLL <- 0)
    rTemp = try(uniroot(O1, c(RR0 - eps, 1e9)), silent = T)
    ifelse (!inherits(rTemp, "try-error"), rUL <- rTemp$root, rUL <- Inf)
  } else if (y1 == n1 & y2 == n2) {
    rTemp <- try(uniroot(O1, c(eps, RR0 + eps)), silent = T)
    ifelse (!inherits(rTemp, "try-error"), rLL <- rTemp$root, rLL <- 0)
    rTemp <- try(uniroot(O1, c(RR0 - eps, min(1 + 1/rLL, 1e+09))), silent = T)
    ifelse (!inherits(rTemp, "try-error"), rUL <- rTemp$root, rUL <- Inf)
  }

  CI.se <- sqrt(1/y1 - 1/n1 + 1/y2 - 1/n2)
  CI.low <- RR0 / exp(qnorm(0.5 + conf.level/2) * CI.se)
  CI.up <- RR0 * exp(qnorm(0.5 + conf.level/2) * CI.se)

  demo <- c("Total Counts" = n, "Conf. Level" = conf.level,
            "Cutoff Value k" = exp(logk), "maxLL" = maxLL)
  contin_2x2 <- matrix(c(y1, y2, y1 + y2,
                         d1, d2, d1 + d2,
                         n1, n2, n,
                         p1, p2, RR0), nrow = 3)
  colnames(contin_2x2) = c("Event (+)", "Event (-)", "Total", "Proportion")
  rownames(contin_2x2) = c("Exposure (+)", "Exposure (-)",
                           "Total / Relative Risk")
  LI <- c("Point Estimate" = RR0, "lower" = rLL, "upper" = rUL, "width" = rUL - rLL)
  CI <- c("Point Estimate" = RR0, "lower" = CI.low, "upper" = CI.up, "width" = CI.up - CI.low)

  plot <- c("Current plot setting is" = plot)
  z <- list(demo = demo, contin_2x2 = contin_2x2, LI = LI, CI = CI, plot = plot)

  plot2 <- c("! Plot CANNOT be drawn :: DATA not supported")
  z2 <- list(demo = demo, contin_2x2 = contin_2x2, LI = LI, CI = CI, plot = plot2)

  if (y1 < n1 & y2 < n2) {
    drr <- seq(0, rUL * 2, length.out = 1e3)
    if (plot == "all") {
      O1plot <- function(rr) {
        plot(rr, O1t(rr), type = "l",
             xlab = "Relative Risk Value",
             ylab = "maxLL - LL - logk",
             main = "Adj. Objective Function (O1 type)")
        abline(h = 0, col = "red")
        abline(v = RR0, lty=2)
        legend("topright",
               legend = c(paste("PE = ", format(RR0, digits = 2)), "Zero Line"),
               lty = c(2, 1),
               col = c("black", "red"))
      }
      O2plot <- function(rr) {
        plot(rr, O2t(rr), type = "l",
             xlab = "Relative Risk Value",
             ylab = "LL",
             main = "Adj. Log Likelihood Function (O2 type)")
        abline(h = maxLL, col = "blue")
        abline(v = RR0, lty=2)
        abline(h = maxLL - logk, col = "red")
        legend("bottomright",
               legend = c(paste("PE = ", format(RR0, digits=2)),
                          paste("maxLL = ", format(maxLL, digits=4)),
                          paste("maxLL-logk = ", format(maxLL-logk, digits=4))),
               lty = c(2, 1, 1),
               col = c("black", "blue", "red"))
      }
      O3plot <- function(rr) {
        plot(rr, O3t(rr), type = "l",
             xlab = "Relative Risk Value",
             ylab = "maxLL - LL",
             main = "Adj. Log LRT (O3 type)")
        abline(h = logk, col = "red")
        abline(v = RR0, lty = 2)
        legend("topright",
               legend = c(paste("PE = ", format(RR0, digits = 2)),
                          paste("logk = ", format(logk, digits = 4))),
               lty = c(2, 1),
               col = c("black", "red"))
      }
      par(mfrow = c(2,2))
      O1plot(drr); O2plot(drr); O3plot(drr)
      par(mfrow = c(1,1))
    } else if (plot == "OBJ" | plot == "O1" | plot == 1) {
      O1plot(drr)
    } else if (plot == "OFV" | plot == "O2" | plot == 2) {
      O2plot(drr)
    } else if (plot == "LRT" | plot == "O3" | plot == 3) {
      O3plot(drr)
    } else {}
    z
  }
  else return(z2)
}
