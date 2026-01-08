#' @importFrom stats binom.test nlminb poisson.test qchisq qf qnorm qt t.test uniroot var var.test
#' @importFrom graphics abline legend par points
#' @export

liod = function (exposure.O_event.O, exposure.O_TOTAL,
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
  o1 <- y1/d1
  o2 <- y2/d2

  n <- n1 + n2
  OR0 <- y1/y2 * d2/d1

  ifelse (y1 == 0 | y1 == n1 | y2 == 0 | y2 == n2, maxLL <- 0,
          maxLL <- y1 * log(p1) + d1 * log(1 - p1) + y2 * log(p2) + d2 * log(1 - p2))


  if (!missing(k)) {
    logk <- log(k)
  } else if (n == 1) {
    logk <- log(2/(1 - conf.level))
  } else {
    logk <- n/2 * log(1 + qf(conf.level, 1, n - 1)/(n - 1))
    logk <- min(logk, log(2/(1 - conf.level)))
  }

  O1 <- function(or) {
    A <- n2 * (or - 1)
    B <- n1 * or + n2 - (y1 + y2) * (or - 1)
    C1 <- -(y1 + y2)
    p2t <- (-B + sqrt(B * B - 4 * A * C1))/(2 * A)
    p1t <- p2t * or/(1 + p2t * (or - 1))
    ifelse (p1t < eps | p1t > 1 - eps | p2t < eps | p2t > 1 - eps, LL <- 0,
            LL <- y1 * log(p1t) + d1 * log(1 - p1t) + y2 * log(p2t) + d2 * log(1 - p2t))
    maxLL - LL - logk
  }

  LLt <- function(or) {
    A <- n2 * (or - 1)
    B <- n1 * or + n2 - (y1 + y2) * (or - 1)
    C1 <- -(y1 + y2)
    p2t <- (-B + sqrt(B * B - 4 * A * C1))/(2 * A)
    p1t <- p2t * or/(1 + p2t * (or - 1))
    y1 * log(p1t) + d1 * log(1 - p1t) + y2 * log(p2t) + d2 * log(1 - p2t)
  }

  O1t <- function(or) maxLL - LLt(or) - logk
  O2t <- function(or) LLt(or)
  O3t <- function(or) maxLL - LLt(or)

  if (y1 == 0 & y2 == 0) {
    oLL <- 0
    oUL <- Inf
  } else if (y1 == 0 & y2 < n2) {
    oLL <- 0
    rTemp <- try(uniroot(O1, c(eps, 1e+09)), silent = T)
    ifelse (!inherits(rTemp, "try-error"), oUL <- rTemp$root, oUL <- Inf)
  } else if (y1 == 0 & y2 == n2) {
    oLL <- 0
    rTemp <- try(uniroot(O1, c(eps, 1e+09)), silent = T)
    ifelse (!inherits(rTemp, "try-error"), oUL <- rTemp$root, oUL <- Inf)
  } else if (y1 < n1 & y2 == 0) {
    rTemp <- try(uniroot(O1, c(eps, 1e+09)), silent = T)
    ifelse (!inherits(rTemp, "try-error"), oLL <- rTemp$root, oLL <- 0)
    oUL <- Inf
  } else if (y1 < n1 & y2 < n2) {
    rTemp <- try(uniroot(O1, c(1e-04, OR0)), silent = T)
    ifelse (!inherits(rTemp, "try-error"), oLL <- rTemp$root, oLL <- 0)
    rTemp <- try(uniroot(O1, c(OR0, 10000)), silent = T)
    ifelse (!inherits(rTemp, "try-error"), oUL <- rTemp$root, oUL <- Inf)
  } else if (y1 < n1 & y2 == n2) {
    oLL <- 0
    rTemp <- try(uniroot(O1, c(eps, 1e+09)), silent = T)
    ifelse (!inherits(rTemp, "try-error"), oUL <- rTemp$root, oUL <- Inf)
  } else if (y1 == n1 & y2 == 0) {
    rTemp <- try(uniroot(O1, c(eps, 1e+09)), silent = T)
    ifelse (!inherits(rTemp, "try-error"), oLL <- rTemp$root, oLL <- 0)
    oUL <- Inf
  } else if (y1 == n1 & y2 < n2) {
    rTemp <- try(uniroot(O1, c(eps, 1e+09)), silent = T)
    ifelse (!inherits(rTemp, "try-error"), oLL <- rTemp$root, oLL <- 0)
    oUL <- Inf
  } else if (y1 == n1 & y2 == n2) {
    oLL <- 0
    oUL <- Inf
  }

  CI.se <- sqrt(1/y1 + 1/y2 + 1/d1 + 1/d2)
  CI.low <- OR0 / exp(qnorm(0.5 + conf.level/2) * CI.se)
  CI.up <- OR0 * exp(qnorm(0.5 + conf.level/2) * CI.se)
  demo <- c("Total Counts" = n, "Conf. Level" = conf.level,
            "Cutoff Value k" = exp(logk), "maxLL" = maxLL)
  contin_2x2 <- matrix(c(y1, y2, y1 + y2,
                         d1, d2, d1 + d2,
                         n1, n2, n,
                         y1/d1, y2/d2, OR0), nrow = 3)
  colnames(contin_2x2) = c("Event (+)", "Event (-)", "Total", "Odds")
  rownames(contin_2x2) = c("Exposure (+)", "Exposure (-)", "Total / Odds Ratio")
  LI <- c("Point Estimate" = OR0, "lower" = oLL, "upper" = oUL, "width" = oUL - oLL)
  CI <- c("Point Estimate" = OR0, "lower" = CI.low, "upper" = CI.up, "width" = CI.up - CI.low)

  plot <- c("Current plot setting is" = plot)
  z <- list(demo = demo, contin_2x2 = contin_2x2, LI = LI, CI = CI, plot = plot)

  plot2 <- c("! Plot CANNOT be drawn :: DATA not supported")
  z2 <- list(demo = demo, contin_2x2 = contin_2x2, LI = LI, CI = CI, plot = plot2)

  if (y1 < n1 & y2 < n2) {
    dor <- seq(0, 2 * oUL, length.out = 1e5)
    if (plot == "all") {
      O1plot <- function(or) {
        plot(or, O1t(or), type = "l",
             xlab = "Odds Ratio Value",
             ylab = "maxLL - LL - logk",
             main = "Adj. Objective Function (O1 type)")
        abline(h = 0, col = "red")
        abline(v = OR0, lty=2)
        legend("topright",
               legend = c(paste("PE = ", format(OR0, digits = 2)), "Zero Line"),
               lty = c(2, 1),
               col = c("black", "red"))
      }
      O2plot <- function(or) {
        plot(or, O2t(or), type = "l",
             xlab = "Odds Ratio Value",
             ylab = "LL",
             main = "Adj. Log Likelihood Function (O2 type)")
        abline(h = maxLL, col = "blue")
        abline(v = OR0, lty=2)
        abline(h = maxLL - logk, col = "red")
        legend("bottomright",
               legend = c(paste("PE = ", format(OR0, digits=2)),
                          paste("maxLL = ", format(maxLL, digits=4)),
                          paste("maxLL-logk = ", format(maxLL-logk, digits=4))),
               lty = c(2, 1, 1),
               col = c("black", "blue", "red"))
      }
      O3plot <- function(or) {
        plot(or, O3t(or), type = "l",
             xlab = "Odds Ratio Value",
             ylab = "maxLL - LL",
             main = "Adj. Log LRT (O3 type)")
        abline(h = logk, col = "red")
        abline(v = OR0, lty = 2)
        legend("topright",
               legend = c(paste("PE = ", format(OR0, digits = 2)),
                          paste("logk = ", format(logk, digits = 4))),
               lty = c(2, 1),
               col = c("black", "red"))
      }
      par(mfrow = c(2,2))
      O1plot(dor); O2plot(dor); O3plot(dor)
      par(mfrow = c(1,1))
    } else if (plot == "OBJ" | plot == "O1" | plot == 1) {
      O1plot(dor)
    } else if (plot == "OFV" | plot == "O2" | plot == 2) {
      O2plot(dor)
    } else if (plot == "LRT" | plot == "O3" | plot == 3) {
      O3plot(dor)
    } else {}
    z
  }
  else return(z2)
}
