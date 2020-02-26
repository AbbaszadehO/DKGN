//' Get gartial correlation coefiicients from given precision matrix
//' @name PCorr
//' @param precision Precision matrix.
//'
//' @return A partial correlation matrix
//' @export

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
NumericMatrix pCorr(NumericMatrix precision) {
  size_t p = precision.nrow();
  NumericMatrix pCor(p,p);
  NumericMatrix p_precision(p,p);

  for(size_t i = 0; i < p-1; i++){
    for(size_t j =i+1; j < p; j++){
      pCor(i,j) = -precision(i,j) / sqrt(precision(i,i)*precision(j,j));
      pCor(j,i) = pCor(i,j);
    }
  }

  return pCor;
}
