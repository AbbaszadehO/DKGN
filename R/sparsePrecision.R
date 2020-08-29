#' L1 regularization to estimate sparse precision matrix
#'
#' @param S Covariance Matrix
#' @param n Number of the samples
#' @param numTF Number of the Transcription Factors
#' @param max_iter Maximum number of iteration in the ADMM algorithm
#' @param rho Penalty parameter of the augmented Lagrangian function
#' @param tol The tolerance for the optimization
#'
#' @return Sparse precision matrix
#' @export
#' @examples

sparsePrecision <-
  function(S,
           n,
           numTF,
           gamma,
           rho,
           max_iter,
           tol
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
    if (numTF > dim(S)[1]) {
      stop("DKGN : Number of trancription factor must be lower than the number of genes. ")
    }
    if (!is.numeric(numTF)) {
      stop("DKGN : numTF must be numeric. ")
    }
    if (!is.numeric(gamma)) {
      stop("DKGN : gamma should be numeric. ")
    }
    if (dim(gamma)[1] != dim(gamma)[2]) {
      stop("DKGN : gamma matrix must be square matrix. ")
    }
    Chat = ADMM(S, gamma, rho, max_iter = max_iter, tol, numvariable = numTF)
    precision = matrix(Chat, nrow = nrow(S))
    return(precision)
  }
