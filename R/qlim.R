#' @importFrom stats binom.test nlminb poisson.test qchisq qf qnorm qt t.test uniroot var var.test
#' @importFrom graphics abline legend par points
#' @export

qlim = function (set.a, set.b, conf.level = 0.95, df = 1, k)
{
  x <- set.a[!is.na(set.a)]
  y <- set.b[!is.na(set.b)]
  n1 <- length(x)
  n2 <- length(y)
  n0 <- n1 + n2

  if (!is.numeric(x) | sum(is.infinite(x) > 0) | sum(is.nan(x)) >
      0 | n1 < 3 | length(unique(x)) == 1)
    stop("Check the input!")
  if (!is.numeric(y) | sum(is.infinite(y) > 0) | sum(is.nan(y)) >
      0 | n2 < 3 | length(unique(y)) == 1)
    stop("Check the input!")

  m1 <- mean(x)
  m2 <- mean(y)
  m0 <- m1 - m2

  if (!missing(k)) {
    logk <- log(k)
  } else {
    logk <- n0/2 * log(1 + qf(conf.level, 1, n0 - df)/(n0 - df))
    logk <- min(logk, log(2/(1 - conf.level)))
  }

  LL1 <- function(th) -n1/2 * (log(2 * pi * sum((x - th)^2)/n1) + 1)
  LL2 <- function(th) -n2/2 * (log(2 * pi * sum((y - th)^2)/n2) + 1)
  maxLL1 <- sapply(m1, LL1)
  maxLL2 <- sapply(m2, LL2)
  maxLL <- sapply(m1, LL1) + sapply(m2, LL2)

  O2 <- function(d) {
    mLLt <- function(m2t) -(sapply(m2t + d, LL1) + sapply(m2t, LL2))
    m2t <- nlminb(start = m2, objective = mLLt)
    -m2t$objective
  }
  O1 <- function(d) maxLL - O2(d) - logk

  meanLL <- uniroot(O1, c(-1e4, m0))$root
  meanUL <- uniroot(O1, c(m0, 1e8))$root

  CI <- t.test(x, y, conf.level = conf.level)
  demo <- c("size a" = n1, "mean a" = m1, "size b" = n2, "mean b" = m2,
            "Conf. Level" = conf.level, "LI df" = df)
  interval <- c("LI PE (= CI PE)" = m0,
                "LI lower" = meanLL, "LI upper" = meanUL,
                "LI width" = meanUL-meanLL,
                "CI lower" = CI$conf.int[1], "CI upper" = CI$conf.int[2],
                "CI width" = CI$conf.int[2]-CI$conf.int[1])
  stat <- c("p-value" = CI$p.value, "Cutoff Value k" = exp(logk),
            "Likelihood Ratio" = exp(O2(0))/exp(maxLL),
            "HETEROsced" = ifelse (O1(0) < 0, FALSE, TRUE))
  verdict <- ifelse (O1(0) < 0, "Equal Mean :: by LI method",
                     "! INequivalent Mean :: by LI method")

  if (0 < meanLL) {
    dm <- seq(-1, meanUL + 1, length.out = 1e3)
  } else if (meanUL < 0) {
    dm <- seq(meanLL - 1, 1, length.out = 1e3)
  } else {
    dm <- seq(meanLL - 1, meanUL + 1, length.out = 1e3)
  }

  plot(dm, sapply(dm, O2), type = "l",
       xlab = "Mean Difference Value",
       ylab = "LL",
       main = "Log Likelihood Function (O2 type)")
  points(0, O2(0), col = "green4", cex = 1.5, pch = 16)
  abline(h = maxLL, col = "blue")
  abline(v = m0, lty = 3)
  abline(h = maxLL - logk, col = "red")
  abline(v = 0, col = "green4", lty = 2)
  abline(h = O2(0), col = "green4")
  legend("bottom", horiz = TRUE, cex = 0.6,
         legend = c(paste("PE = ", format(m0, digits = 2)),
                    paste("maxLL = ", format(maxLL, digits = 4)),
                    paste("Cutoff = ", format(maxLL-logk, digits = 4)),
                    paste("LL Eq Mean = ", format(O2(0), digits = 4))),
         lty = c(3, 1, 2, 1),
         col = c("black", "blue", "red", "green4"))

  return(list(demo = demo, interval = interval, stat = stat, verdict = verdict))
}
