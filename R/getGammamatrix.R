#' Gamma matrix for L1 regularization
#'
#' @param S Covariance Matrix
#' @param confidence Confidence level
#' @param prior Prior knowledge matrix
#'
#' @return Gamma matrix
#' @export
#' @examples

getGammamatrix <-
  function(S,
           confidence,
           prior=NULL
  ) {
    if (dim(S)[1] != dim(S)[2]) {
      stop("DKGN : S must be square matrix. ")
    }
    if (!is.numeric(S)) {
      stop("DKGN : S must be numeric. ")
    }
    if (any(is.na(S)) || any(is.nan(S))) {
      stop("DKGN: S contains missing values. ")
    }
    if (!is.numeric(confidence)) {
      stop("DKGN : confidence should be numeric. ")
    }
    if(!is.null(prior)){
      if (dim(prior)[1] != dim(prior)[2]) {
        stop("DKGN : Prior matrix must be square matrix. ")
      }
      if (!is.numeric(prior)) {
        stop("DKGN : Prior must be numeric. ")
      }
      if (any(is.na(prior)) || any(is.nan(prior))) {
        stop("DKGN: Prior matrix contains missing values. ")
      }
    }
    if(is.null(prior)){
      prior = matrix(data = 1, nrow = nrow(S), ncol = ncol(S))
    }
    diagS = (diag(S))
    SS = outer(diagS, diagS)
    maxSS = max(SS)
    tval  = qt(0.5 + (confidence / 2), df = (n - 2))
    gammav   = tval * maxSS / sqrt(n - 2 + (tval ^ 2))
    gamma = matrix(gammav,nrow = nrow(S),ncol = ncol(S))
    if(!is.null(prior)){
      for (i in seq(nrow(S))){
        for (j in seq(nrow(S))){
          gamma[i,j] =   prior[i,j] * tval * maxSS / sqrt(n - 2 + (tval ^ 2))
        }
      }
    }
    return(gamma)
  }
