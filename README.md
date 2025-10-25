
# inteli

<!-- badges: start -->
<!-- badges: end -->

Provide likelihood-based functions and interval estimation results, compared to
the traditional statistics (CI). This is the expanded version for {LBI}- and 
{wnl}-package, formulated by Kyun-Seop Bae <k@acr.kr>.


## Installation

You can install the development version of inteli from [GitHub](https://github.com/) with:

``` r
# install.packages("pak")
pak::pak("minQk-acr/inteli")
```

## Example

This is a basic example which shows you how to solve a common problem:

``` r
library(inteli)

x <- rnorm(20, 0, 1)
y <- rnorm(40, 0, 1)

varE(x)
varE(y)
varR(x, y)

varEplot(x, FALSE)
varEplot(x, TRUE)

varEplot(y, FALSE)
varEplot(y, TRUE)

varRplot(x, y, FALSE)
varRplot(x, y, TRUE)
```

