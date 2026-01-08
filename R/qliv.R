#' @importFrom stats binom.test nlminb poisson.test qchisq qf qnorm qt t.test uniroot var var.test
#' @importFrom graphics abline legend par points
#' @export

qliv = function (set.a, set.b, conf.level = 0.95, df = 2.4, k)
{
  dx <- set.a[!is.na(set.a)]
  dy <- set.b[!is.na(set.b)]
  dn1 <- length(dx)
  dn2 <- length(dy)

  if (dn1 <= dn2) {
    x <- dx
    y <- dy
    } else {
      y <- dx
      x <- dy
    }

  n1 <- length(x)
  n2 <- length(y)
  n1v1 <- sum((x - mean(x))^2)
  n2v2 <- sum((y - mean(y))^2)
  v1 <- n1v1/n1
  v2 <- n2v2/n2

  if (!is.numeric(x) | sum(is.infinite(x) > 0) | sum(is.nan(x)) >
      0 | n1 < 3 | length(unique(x)) == 1)
    stop("Check the input!")
  if (!is.numeric(y) | sum(is.infinite(y) > 0) | sum(is.nan(y)) >
      0 | n2 < 3 | length(unique(y)) == 1)
    stop("Check the input!")

  R0 <- v1/v2
  n0 <- n1 + n2
  maxLL <- -(n1 * (log(2 * pi * v1) + 1) + n2 * (log(2 * pi * v2) + 1))/2

  if (!missing(k)) {
    logk <- log(k)
    } else {
      logk <- n0/2 * log(1 + qf(conf.level, 1, n0 - df)/(n0 - df))
      logk <- min(logk, log(2/(1 - conf.level)))
    }

  O2 <- function(r) {
    th <- (n1v1 + r * n2v2) / r / n0
    ln2pith <- log(2 * pi * th)
    -(n1 * (log(r) + ln2pith) + n1v1 / r / th + n2 * ln2pith + n2v2 / th) / 2
  }

  O1 <- function(r) maxLL - O2(r) - logk

  varLL <- uniroot(O1, c(1e-08, R0))$root
  varUL <- uniroot(O1, c(R0, 1e+06))$root
  CI <- var.test(x, y, conf.level = conf.level)
  demo <- c("Size 1 (smaller)" = n1, "Size 2" = n2,
            "Conf. Level" = conf.level, "LI df" = df)
  interval <- c("LI PE" = R0,
                "LI lower" = varLL,
                "LI upper" = varUL,
                "LI width" = varUL-varLL,
                "CI PE" = unname(CI$estimate),
                "CI lower" = CI$conf.int[1],
                "CI upper" = CI$conf.int[2],
                "CI width" = CI$conf.int[2]-CI$conf.int[1])
  stat <- c("p-value" = CI$p.value,
            "Cutoff Value k" = exp(logk),
            "Likelihood Ratio" = exp(O2(1))/exp(maxLL),
            "HETEROscedasticity" = ifelse (O1(1) < 0, FALSE, TRUE))
  verdict <- ifelse (O1(1) < 0,
                     "Equal Variance :: by LI method",
                     "! INequivalent Variance :: by LI method")


  if (varUL * 2 < 1) {
    dr <- seq(varLL/2, 2, length.out = 1e3)
  } else if (1 < varLL/2) {
    dr <- seq(0.5, varUL * 2, length.out = 1e3)
  } else {dr <- seq(varLL/2, varUL * 2, length.out = 1e3)}

  plot(dr, O2(dr), type = "l",
       xlab = "Variance Ratio Value (Sample Size Allocated)",
       ylab = "LL",
       main = "Log Likelihood Function (O2 type)")
  points(1, O2(1), col = "green4", cex = 1.5, pch = 16)
  abline(h = maxLL, col = "blue")
  abline(v = R0, lty = 3)
  abline(h = maxLL - logk, col = "red")
  abline(v = 1, col = "green4", lty = 2)
  abline(h = O2(1), col = "green4")
  legend("bottom", horiz = TRUE, cex = 0.6,
         legend = c(paste("PE = ", format(R0, digits = 2)),
                    paste("maxLL = ", format(maxLL, digits = 4)),
                    paste("Cutoff = ", format(maxLL-logk, digits = 4)),
                    paste("LL Eq Var = ", format(O2(1), digits = 4))),
         lty = c(3, 1, 2, 1),
         col = c("black", "blue", "red", "green4"))

  return(list(demo = demo, interval = interval, stat = stat, verdict = verdict))
}
