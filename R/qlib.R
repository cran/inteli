#' @importFrom stats binom.test nlminb poisson.test qchisq qf qnorm qt t.test uniroot var var.test
#' @importFrom graphics abline legend par points
#' @export

qlib = function (event, total, test.val, conf.level = 0.95, eps = 1e-08, k)
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

  if (p0 < eps) {
    bLL <- 0
  } else {
    rTemp <- try(uniroot(O1, c(max(eps, p0 - LE), p0 + eps)), silent = TRUE)
    ifelse (!inherits(rTemp, "try-error"), bLL <- rTemp$root, bLL <- max(p0 - LE, 0))
  }

  if (p0 > 1 - eps) {
    bUL <- 1
  } else {
    rTemp = try(uniroot(O1, c(p0 - eps, min(p0 + LE, 1 - eps))), silent = TRUE)
    ifelse (!inherits(rTemp, "try-error"), bUL <- rTemp$root, bUL <- min(p0 + LE, 1))
  }

  CI <- binom.test(event, total, test.val, conf.level = conf.level)
  demo <- c("Event Counts" = y, "Total Counts" = n,
            "Conf. Level" = conf.level, "Testing Value" = test.val)
  interval <- c("LI PE (= CI PE)" = p0,
                "LI lower" = bLL,
                "LI upper" = bUL,
                "LI width" = bUL-bLL,
                "CI lower" = CI$conf.int[1],
                "CI upper" = CI$conf.int[2],
                "CI width" = CI$conf.int[2]-CI$conf.int[1])
  stat <- c("p-value" = CI$p.value,
            "Cutoff Value k" = exp(logk),
            "Likelihood Ratio" = exp(O2(test.val))/exp(maxLL),
            "Testing Result" = ifelse(O1(test.val) < 0, TRUE, FALSE))
  verdict <- ifelse (O1(test.val) < 0,
                     "Testing value can be explained :: by LI method",
                     "! Testing value !CANNOT! be explained :: by LI method")

  dp <- seq(0, 1, length.out = 1e3)
  plot(dp, O2(dp), type = "l",
       xlab = "Proportion Value",
       ylab = "LL",
       main = "Log Likelihood Function (O2 type)")
  points(test.val, O2(test.val), col = "green4", cex = 1.5, pch = 16)
  abline(h = maxLL, col = "blue")
  abline(v = p0, lty = 3)
  abline(h = maxLL - logk, col = "red")
  abline(v = test.val, col = "green4", lty = 2)
  abline(h = O2(test.val), col = "green4")
  legend("bottomright",
         legend = c(paste("PE = ", format(p0, digits = 2)),
                    paste("maxLL = ", format(maxLL, digits = 4)),
                    paste("Cutoff = ", format(maxLL-logk, digits = 4)),
                    paste("Testing = ", format(test.val, digits = 4)),
                    paste("LL Testing = ", format(O2(test.val), digits = 4))),
         lty = c(3, 1, 1, 2, 1),
         col = c("black", "blue", "red", "green4", "green4"))

  return(list(demo = demo, interval = interval, stat = stat, verdict = verdict))
}





