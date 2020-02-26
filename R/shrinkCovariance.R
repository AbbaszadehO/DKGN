#' Linear shrinkage of sample covariance matrix and a target matrix based on empirical Bayesian approach
#'
#' @param S Sample covariance matrix
#' @param target Target matrix
#' @param n Number of the samples
#' @param lambda A sequence between [0,1]
#'
#' @return A shrunken covariance matrix
#' @export
#'
#' @examples

shrinkCovariance <-
  function(S, target, n, lambda = seq(0.01, 0.99, 0.01)) {
    if (dim(S)[1] != dim(S)[2]) {
      stop("DKGN : S must be square matrix. ")
    }
    if (!is.numeric(S)) {
      stop("DKGN : S must be numeric. ")
    }
    if (any(is.na(S)) || any(is.nan(S))) {
      stop("DKGN: S contains missing values. ")
    }

    if (dim(target)[1] != dim(target)[2]) {
      stop("DKGN : target must be square matrix. ")
    }
    if (!is.numeric(target)) {
      stop("DKGN : target must be numeric. ")
    }
    if (any(is.na(target)) || any(is.nan(target))) {
      stop("DKGN: target contains missing values. ")
    }

    lambda_star = getIntensity(S, target = target, lambda = lambda, n)
    message(paste("DKGN: optimal shrinkage intensity is: ", lambda_star, sep = ""))
    sigma_hat = lambda_star * target + (1 - lambda_star) * S
    return(sigma_hat)
  }
