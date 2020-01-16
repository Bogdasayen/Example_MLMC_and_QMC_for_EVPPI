##
##  UCS / The R ToolBox
##

## Filename: sfunc.R
## Modified: Sat Feb 21 18:24:49 2004 (evert)   
##   Author: Stefan Evert
##  Purpose: Special functions (beta, gamma, confidence intervals)


## (complete) gamma function and its logarithm (all logarithms are base 10)
Cgamma <- function (a, log=FALSE) {
  if (log) {
    lgamma(a) / log(10)
  } else {
    gamma(a)
  }
}

## regularised gamma function (lower P(a,x) and upper Q(a,x)) and its inverse
Rgamma <- function (a, x, lower=TRUE, log=FALSE) {
  if (log) {
    pgamma(x, shape=a, scale=1, lower.tail=lower, log=TRUE) / log(10)
  } else {
    pgamma(x, shape=a, scale=1, lower.tail=lower, log=FALSE)
  }
}
Rgamma.inv <- function(a, y, lower=TRUE, log=FALSE) {
  if (log) {
    qgamma(y * log(10), shape=a, scale=1, lower.tail=lower, log=TRUE)
  } else {
    qgamma(y, shape=a, scale=1, lower.tail=lower, log=FALSE)
  }
}

## incomplete gamma function (lower gamma(a,x) and upper Gamma(a,x)) and its inverse
Igamma <- function (a, x, lower=TRUE, log=FALSE) {
  if (log) {
    Cgamma(a, log=TRUE) + Rgamma(a, x, lower, log=TRUE)
  } else {
    Cgamma(a, log=FALSE) * Rgamma(a, x, lower, log=FALSE)
  }
}
Igamma.inv <- function (a, y, lower=TRUE, log=FALSE) {
  if (log) {
    Rgamma.inv(a, y - Cgamma(a, log=TRUE), lower, log=TRUE)
  } else {
    Rgamma.inv(a, y / Cgamma(a, log=FALSE), lower, log=FALSE)
  }
}

## beta function and its logarithm
Cbeta <- function(a, b, log=FALSE) {
  if (log) {
    lbeta(a, b) / log(10)
  } else {
    beta(a, b)
  }
}

## regularised beta function I(x; a, b) and its inverse
Rbeta <- function (x, a, b, log=FALSE) {
  if (log) {
    pbeta(x, shape1=a, shape2=b, log=TRUE) / log(10)
  } else {
    pbeta(x, shape1=a, shape2=b, log=FALSE)
  }
}
Rbeta.inv <- function (y, a, b, log=FALSE) {
  if (log) {
    qbeta(y * log(10), shape1=a, shape2=b, log=TRUE)
  } else {
    qbeta(y, shape1=a, shape2=b, log=FALSE)
  }
}

## incomplete beta function B(x; a, b) and its inverse
Ibeta <- function (x, a, b, log=FALSE) {
  if (log) {
    Cbeta(a, b, log=TRUE) + Rbeta(x, a, b, log=TRUE)
  } else {
    Cbeta(a, b, log=FALSE) * Rbeta(x, a, b, log=FALSE)
  }
}
Ibeta.inv <- function (y, a, b, log=FALSE) {
  if (log) {
    Rbeta.inv(y - Cbeta(a, b, log=TRUE), a, b, log=TRUE)
  } else {
    Rbeta.inv(y / Cbeta(a, b, log=FALSE), a, b, log=FALSE)
  }
}

## two-sided confidence interval for success probability p of binomial distribution,
## given that k successes out of size trials have been observed
##   p.limit <- binom.conf.interval(k, size, conf.level=0.05, limit="lower", one.sided=FALSE)
## conf.level = confidence level, e.g. 0.01 corresponds to 99% confidence
## one.sided  = if TRUE, one-sided confidence level, otherwise two-sided
binom.conf.interval <- function(k, size, limit=c("lower","upper"), conf.level=0.05, one.sided=FALSE) {
  l <- length(k)
  if (l != length(size)) {
    if (length(size)==1) {
      size <- rep(size, times=l)
    } else {
      stop("Parameters k and size must be vectors of identical length.")
    }
  }
  limit <- match.arg(limit)
  # use regularised incomplete Beta function (= distribution function of Beta distribution)
  # to compute two-sided confidence intervals for parameter of binomial distribution
  if (one.sided) alpha <- conf.level else alpha <- conf.level / 2
  if (limit == "lower") {
    return(qbeta(alpha, k, size - k + 1))
  }
  else {
    return(qbeta(1 - alpha, k + 1, size - k))
  }
}
