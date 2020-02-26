#' Get target matrix for Covariance Estimation
#'
#' @param S Sample Covariance Matrix
#' @param type  Dteremine the target matrix type. default = "mean"
#' @param prior A Prior matrix with dimensions equal to S.
#' @return A target matrix
#' @export
#'
#' @examples
#'
getTarget <- function(S,
                      type = c("mean", "identity", "variance", "prior"),
                      prior = NULL) {
  if (dim(S)[1] != dim(S)[2]) {
    stop("KDGN : S must be square matrix. ")
  }
  if (!is.numeric(S)) {
    stop("KDGN : S must be numeric. ")
  }
  if (any(is.na(S)) || any(is.nan(S))) {
    stop("KDGN: S contains missing values. ")
  }

  sh_type = match.arg(type)


  p = nrow(S)
  target = matrix(0, nrow = p, ncol = p)

  if (sh_type == "identity") {
    target = diag(1, p)
  }

  if (sh_type == "variance") {
    target = diag(diag(S), p)
  }

  if (sh_type == "mean") {
    target = diag(mean(diag(S)), p)
  }

  if (sh_type == "prior") {
    prior_one = diag(1, p)
    rownames(prior_one) = rownames(S)
    colnames(prior_one) = rownames(S)
    for (i in 1:dim(prior)[1]) {
      node1 = prior[i, 1]
      if (!is.element(node1, row.names(S))) {
        stop("KDGN: Gene name not found. ")
      }
      node2 = prior[i, 2]
      if (!is.element(node2, row.names(S))) {
        stop("KDGN: Gene name not found. ")
      }
      prior_one[node1, node2] = 1
      prior_one[node2, node1] = 1
    }
    for (i in 1:p) {
      for (j in 1:i) {
        if (i != j &
            prior_one[i, j] == 1)
          target[i, j] = S[i, j]
        if (i == j)
          target[i, j] = S[i, j]
        target[j, i] = target[i, j]
      }
    }
  }
  return(target)
}
