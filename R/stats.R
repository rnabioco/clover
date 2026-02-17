# Statistical utility functions -------------------------------------------------

#' Calculate log2 fold change with pseudocount.
#'
#' @param x Numeric vector of numerator values.
#' @param y Numeric vector of denominator values.
#' @param pseudocount Pseudocount added to both numerator and
#'   denominator before computing the ratio. Default `1`.
#'
#' @return Numeric vector of log2 fold changes.
#'
#' @export
#'
#' @examples
#' calc_fold_change(10, 5)
#' calc_fold_change(c(0, 10), c(10, 0), pseudocount = 1)
calc_fold_change <- function(x, y, pseudocount = 1) {
  log2((x + pseudocount) / (y + pseudocount))
}

#' Calculate Cohen's d effect size.
#'
#' Compute the standardized mean difference between two groups using
#' the pooled standard deviation.
#'
#' @param x Numeric vector for the first group.
#' @param y Numeric vector for the second group.
#'
#' @return A single numeric value for the effect size.
#'
#' @export
#'
#' @examples
#' cohens_d(rnorm(20, mean = 5), rnorm(20, mean = 3))
cohens_d <- function(x, y) {
  nx <- length(x)
  ny <- length(y)
  mx <- mean(x, na.rm = TRUE)
  my <- mean(y, na.rm = TRUE)
  sx <- stats::sd(x, na.rm = TRUE)
  sy <- stats::sd(y, na.rm = TRUE)

  pooled_sd <- sqrt(((nx - 1) * sx^2 + (ny - 1) * sy^2) / (nx + ny - 2))
  (mx - my) / pooled_sd
}

#' Propagate standard error for a ratio.
#'
#' Compute the standard error of `a / b` using first-order error
#' propagation.
#'
#' @param a Numeric vector of numerator values.
#' @param b Numeric vector of denominator values.
#' @param se_a Numeric vector of standard errors for `a`.
#' @param se_b Numeric vector of standard errors for `b`.
#'
#' @return Numeric vector of propagated standard errors.
#'
#' @export
#'
#' @examples
#' propagate_error_ratio(10, 5, 1, 0.5)
propagate_error_ratio <- function(a, b, se_a, se_b) {
  abs(a / b) * sqrt((se_a / a)^2 + (se_b / b)^2)
}

#' Propagate standard error for a difference.
#'
#' Compute the standard error of `a - b` assuming independent errors.
#'
#' @param se_a Numeric vector of standard errors for the first value.
#' @param se_b Numeric vector of standard errors for the second value.
#'
#' @return Numeric vector of propagated standard errors.
#'
#' @export
#'
#' @examples
#' propagate_error_diff(1, 0.5)
propagate_error_diff <- function(se_a, se_b) {
  sqrt(se_a^2 + se_b^2)
}
