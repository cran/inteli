#' @importFrom stats binom.test nlminb poisson.test qchisq qf qnorm qt t.test uniroot var var.test
#' @importFrom graphics abline legend par points
#' @export

lird = function (exposure.O_event.O, exposure.O_TOTAL,
                 exposure.X_event.O, exposure.X_TOTAL,
                 plot = "all", conf.level = 0.95, eps = 1e-08, k)
{
  y1 <- exposure.O_event.O
  n1 <- exposure.O_TOTAL
  y2 <- exposure.X_event.O
  n2 <- exposure.X_TOTAL

  if (length(y1) != 1 | length(n1) != 1 |
      (y1 < 0) | (n1 < 0) | !is.finite(y1) | !is.finite(n1))
    stop("Check the input!")
  if (length(y2) != 1 | length(n2) != 1 |
      (y2 < 0) | (n2 < 0) | !is.finite(y2) | !is.finite(n2))
    stop("Check the input!")
  if (n1 > 1/eps | n2 > 1/eps)
    stop("Total size(s) is/are too large")

  p1 <- y1/n1
  p2 <- y2/n2
  d1 <- n1 - y1
  d2 <- n2 - y2
  RD0 <- p1 - p2

  if (p1 == 0) p1 <- eps
  if (p2 == 0) p2 <- eps
  if (p1 == 1) p1 <- 1 - eps
  if (p2 == 1) p2 <- 1 - eps

  maxLL <- y1 * log(p1) + d1 * log(1 - p1) + y2 * log(p2) + d2 * log(1 - p2)
  maxLL <- ifelse(is.finite(maxLL), maxLL, 0)

  n <- n1 + n2

  if (!missing(k)) {
    logk <- log(k)
  } else if (n == 1) {
    logk <- log(2/(1 - conf.level))
  } else {
    logk <- n/2 * log(1 + qf(conf.level, 1, n - 1)/(n - 1))
    logk <- min(logk, log(2/(1 - conf.level)))
  }

  O1 <- function(rd) {
    L3 <- n1 + n2
    L2 <- (n1 + 2 * n2) * rd - L3 - y1 - y2
    L1 <- (n2 * rd - L3 - 2 * y2) * rd + y1 + y2
    L0 <- y2 * rd * (1 - rd)
    q <- L2^3/(3 * L3)^3 - L1 * L2/(6 * L3^2) + L0/(2 * L3)
    if (!is.finite(q) | abs(q) < eps) q <- 0
    p <- sign(q) * sqrt(max(0, L2^2/(3 * L3)^2 - L1/(3 * L3)))
    a <- (pi + ifelse(p == 0, acos(0), acos(min(max(q/p^3, -1), 1))))/3
    p2t <- min(max(eps, 2 * p * cos(a) - L2/(3 * L3)), 1 - eps)
    p1t <- min(max(eps, p2t + rd), 1 - eps)
    LL <- y1 * log(p1t) + d1 * log(1 - p1t) + y2 * log(p2t) + d2 * log(1 - p2t)
    LL <- ifelse(is.finite(LL), LL, 0)
    maxLL - LL - logk
  }

  LLt <- function(rd) {
    L3 <- n1 + n2
    L2 <- (n1 + 2 * n2) * rd - L3 - y1 - y2
    L1 <- (n2 * rd - L3 - 2 * y2) * rd + y1 + y2
    L0 <- y2 * rd * (1 - rd)
    q <- L2^3/(3 * L3)^3 - L1 * L2/(6 * L3^2) + L0/(2 * L3)
    p <- sign(q) * sqrt(max(0, L2^2/(3 * L3)^2 - L1/(3 * L3)))
    a <- (pi + acos(min(max(q/p^3, -1), 1)))/3
    p2t <- 2 * p * cos(a) - L2/(3 * L3)
    p1t <- p2t + rd
    y1 * log(p1t) + d1 * log(1 - p1t) + y2 * log(p2t) + d2 * log(1 - p2t)
  }

  O1t <- function(rd) maxLL - LLt(rd) - logk
  O2t <- function(rd) LLt(rd)
  O3t <- function(rd) maxLL - LLt(rd)

  if (RD0 < -1 + eps) {
    rLL <- -1
  } else {
    rTemp <- try(uniroot(O1, c(-1 + eps, RD0 + eps)), silent = T)
    ifelse (!inherits(rTemp, "try-error"), rLL <- rTemp$root, rLL <- -1)
  }

  if (RD0 > 1 - eps) {
    rUL <- 1
  } else {
    rTemp <- try(uniroot(O1, c(RD0 - eps, 1 - eps)), silent = T)
    ifelse (!inherits(rTemp, "try-error"), rUL <- rTemp$root, rUL <- 1)
  }

  CI.se <- sqrt(p1 * (1 - p1)/n1 + p2 * (1 - p2)/n2)
  CI.low <- RD0 - qnorm(0.5 + conf.level/2) * CI.se
  CI.up <- RD0 + qnorm(0.5 + conf.level/2) * CI.se

  demo <- c("Total Counts" = n, "Conf. Level" = conf.level,
            "Cutoff Value k" = exp(logk), "maxLL" = maxLL)
  contin_2x2 <- matrix(c(y1, y2, y1 + y2,
                         d1, d2, d1 + d2,
                         n1, n2, n,
                         p1, p2, RD0), nrow = 3)
  colnames(contin_2x2) = c("Event (+)", "Event (-)", "Total", "Proportion")
  rownames(contin_2x2) = c("Exposure (+)", "Exposure (-)", "Total / Risk Diff")
  LI <- c("Point Estimate" = RD0, "lower" = rLL, "upper" = rUL, "width" = rUL - rLL)
  CI <- c("Point Estimate" = RD0, "lower" = CI.low, "upper" = CI.up, "width" = CI.up - CI.low)

  plot <- c("Current plot setting is" = plot)
  z <- list(demo = demo, contin_2x2 = contin_2x2, LI = LI, CI = CI, plot = plot)

  plot2 <- c("! Plot CANNOT be drawn :: DATA not supported")
  z2 <- list(demo = demo, contin_2x2 = contin_2x2, LI = LI, CI = CI, plot = plot2)

  if (-1 + eps < RD0 &  RD0 < 1 - eps) {
    drd <- seq(-1, 1, length.out = 1e5)
    if (plot == "all") {
      O1plot <- function(rd) {
        plot(rd, O1t(rd), type = "l",
             xlab = "Risk Difference Value",
             ylab = "maxLL - LL - logk",
             main = "Adj. Objective Function (O1 type)")
        abline(h = 0, col = "red")
        abline(v = RD0, lty=2)
        legend("topright",
               legend = c(paste("PE = ", format(RD0, digits = 2)), "Zero Line"),
               lty = c(2, 1),
               col = c("black", "red"))
      }
      O2plot <- function(rd) {
        plot(rd, O2t(rd), type = "l",
             xlab = "Risk Difference Value",
             ylab = "LL",
             main = "Adj. Log Likelihood Function (O2 type)")
        abline(h = maxLL, col = "blue")
        abline(v = RD0, lty=2)
        abline(h = maxLL - logk, col = "red")
        legend("bottomright",
               legend = c(paste("PE = ", format(RD0, digits=2)),
                          paste("maxLL = ", format(maxLL, digits=4)),
                          paste("maxLL-logk = ", format(maxLL-logk, digits=4))),
               lty = c(2, 1, 1),
               col = c("black", "blue", "red"))
      }
      O3plot <- function(rd) {
        plot(rd, O3t(rd), type = "l",
             xlab = "Risk Difference Value",
             ylab = "maxLL - LL",
             main = "Adj. Log LRT (O3 type)")
        abline(h = logk, col = "red")
        abline(v = RD0, lty = 2)
        legend("topright",
               legend = c(paste("PE = ", format(RD0, digits = 2)),
                          paste("logk = ", format(logk, digits = 4))),
               lty = c(2, 1),
               col = c("black", "red"))
      }
      par(mfrow = c(2,2))
      O1plot(drd); O2plot(drd); O3plot(drd)
      par(mfrow = c(1,1))
    } else if (plot == "OBJ" | plot == "O1" | plot == 1) {
      O1plot(drd)
    } else if (plot == "OFV" | plot == "O2" | plot == 2) {
      O2plot(drd)
    } else if (plot == "LRT" | plot == "O3" | plot == 3) {
      O3plot(drd)
    } else {}
    z
  }
  else return(z2)
}
