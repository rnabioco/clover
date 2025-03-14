# stats -----------------------------------------

#' Kullback-Leibler divergence
#'
#' @examples
#' # test for asymmetry
#' p <- c(0.9, 0.1, 0.0001, 0.0001)
#' q <- c(0.1, 0.9, 0.0001, 0.0001)
#'
#' p_kl <- kl_divergence(p, q)
#' q_kl <- kl_divergence(q, p)
#'
#' abs(p_kl - q_kl)
#'
#' @family stats
#'
#' @export
kl_divergence <- function(p, q) {
  # Avoid log(0) issues
  idx <- which(p > 0)
  sum(p[idx] * log(p[idx] / q[idx]))
}

#' Symmetric Kullback-Leibler divergence
#'
#' @examples
#' p <- c(0.3, 0.2, 0.2, 0.3)
#' q <- c(0.25, 0.25, 0.25, 0.25)
#'
#' symmetric_kl(p, q)
#'
#' @family stats
#'
#' @export
symmetric_kl <- function(p, q) {
  kl_divergence(p, q) + kl_divergence(q, p)
}
