#' L1 regularization to estimate sparse precision matrix
#'
#' @param S Covariance Matrix
#' @param n Number of the samples
#' @param numTF Number of the Transcription Factors
#' @param confidence confidence value
#' @param max_iter Maximum number of iteration in the ADMM algorithm
#' @param rho Penalty parameter of the augmented Lagrangian function
#' @param tol The tolerance for the optimization
#' @param type Type of the ADMM algorithm must be used
#' @param prior Available prior knowledge
#'
#' @return Sparse precision matrix
#' @export
#' @examples

sparsePrecision <-
  function(S,
           n,
           numTF,
           rho,
           confidence,
           max_iter,
           tol,
           type = c("data-driven", "knowledge-based"),
           prior = NULL) {


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
    if (confidence > 1) {
      stop("DKGN : Confidence must be lower than 1. ")
    }
    p_type = match.arg(type)
    if(p_type == "knowledge-based"){
      if (dim(prior)[1] != dim(prior)[2]) {
        stop("DKGN : prior must be square matrix. ")
      }
      if (!is.numeric(prior)) {
        stop("DKGN : prior must be numeric. ")
      }
      if (any(is.na(prior)) || any(is.nan(prior))) {
        stop("DKGN: prior contains missing values. ")
      }
      if(is.null(prior)){
        stop("DKGN : Prior target must be determinded. ")
      }
    }
    diagS = (diag(S))
    SS = outer(diagS, diagS)
    maxSS = max(SS)
    tval  = qt(0.5 + (confidence / 2), df = (n - 2))
    gammav   = tval * maxSS / sqrt(n - 2 + (tval ^ 2))

    if(p_type=="data-driven"){
      Chat = dADMM(S, gammav, rho, max_iter = max_iter, tol, numvariable = numTF)
    }
    else{
      Chat = kADMM(S, gammav, rho, max_iter = max_iter, tol, numvariable = numTF, prior = prior)
    }
    precision = matrix(Chat, nrow = nrow(S))
    return(precision)
  }
