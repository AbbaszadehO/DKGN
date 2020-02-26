//' Function for data standardization
//' @name center
//' @param raw_data A n*p matrix.
//' @param standard A Boolean parameter that determines whether standardization is done.
//'
//' @return A Normal data
//' @export

#include <Rcpp.h>
#include <string>

using namespace Rcpp;
// [[Rcpp::export]]
NumericVector center(NumericMatrix raw_data, bool standard = false) {

  size_t n = raw_data.nrow();
  size_t p = raw_data.ncol();

  NumericMatrix centered_data(n,p);
  for(size_t i = 0; i < p; i++){
    centered_data(_,i) = raw_data(_,i) - mean(raw_data(_,i));
  }
  if (standard == false){
    return centered_data;
  }

  NumericVector l2norm(p);
  NumericMatrix standandard_data(n,p);
  for(size_t i = 0; i < p; i++){
    l2norm[i] = sqrt(sum(pow(centered_data(_,i),2)));
    if(l2norm[i] == 0){
      stop("\n DKGN: Variable " + toString(i+1) + " has zero variance. \n Please remove that.");
    }
    standandard_data(_,i) = centered_data(_,i) / l2norm[i] * sqrt(n);
  }
  return standandard_data;
}
