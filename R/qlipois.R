#' @importFrom stats binom.test nlminb poisson.test qchisq qf qnorm qt t.test uniroot var var.test
#' @importFrom graphics abline legend par points
#' @export

qlipois = function(data, test.val, conf.level = 0.95, eps = 1e-08, k)
{
  x <- data[!is.na(data)]
  n <- length(x)

  if (!is.numeric(x)) stop("! Check the data input!")
  ifelse (n == 1, Lam0 <- x, Lam0 <- sum(x)/n)

  maxLL <- sum(x) * log(Lam0) - n * Lam0

  if (!missing(k)) {
    logk <- log(k)
  } else if (n == 1) {
    logk <- log(2/(1 - conf.level))
  } else {
    logk <- n/2 * log(1 + qf(conf.level, 1, n - 1)/(n - 1))
    logk <- min(logk, log(2/(1 - conf.level)))
  }

  LL <- function(lam, cutoff) ifelse(lam > 0, sum(x) * log(lam) - n * lam, 0)
  O1 <- function(lam, cutoff) maxLL - LL(lam, cutoff) - logk
  O2 <- function(lam, cutoff) LL(lam, cutoff)

  ifelse (Lam0 > 0,
          pLL <- uniroot(O1, c(eps, Lam0), cutoff = logk)$root, pLL <- 0)
  ifelse (!is.finite(maxLL),
          pUL <- Inf, pUL <- uniroot(O1, c(Lam0, 1e+09), cutoff = logk)$root)

  demo <- c("Poisson Mean" = Lam0, "Data Size" = n,
            "Conf. Level" = conf.level, "Testing Value" = test.val)
  interval <- c("LI PE" = Lam0,
                "LI lower" = pLL,
                "LI upper" = pUL,
                "LI width" = pUL-pLL)
  stat <- c("Cutoff Value k" = exp(logk),
            "Likelihood Ratio" = exp(O2(test.val))/exp(maxLL),
            "Testing Result" = ifelse(O1(test.val) < 0, TRUE, FALSE))
  verdict <- ifelse (O1(test.val) < 0,
                     "Testing value can be explained :: by LI method",
                     "! Testing value !CANNOT! be explained :: by LI method")

  dp <- seq(0, 5*Lam0, length.out = 1e3)
  plot(dp, O2(dp), type = "l",
       xlab = "Poisson Mean Value",
       ylab = "LL",
       main = "Log Likelihood Function (O2 type)")
  points(test.val, O2(test.val), col = "green4", cex = 1.5, pch = 16)
  abline(h = maxLL, col = "blue")
  abline(v = Lam0, lty = 3)
  abline(h = maxLL - logk, col = "red")
  abline(v = test.val, col = "green4", lty = 2)
  abline(h = O2(test.val), col = "green4")
  legend("bottomright",
         legend = c(paste("PE = ", format(Lam0, digits = 2)),
                    paste("maxLL = ", format(maxLL, digits = 4)),
                    paste("Cutoff = ", format(maxLL-logk, digits = 4)),
                    paste("Testing = ", format(test.val, digits = 4)),
                    paste("LL Testing = ", format(O2(test.val), digits = 4))),
         lty = c(3, 1, 1, 2, 1),
         col = c("black", "blue", "red", "green4", "green4"))

  return(list(demo = demo, interval = interval, stat = stat, verdict = verdict))
}
