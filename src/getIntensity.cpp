//' Shrinkage estimator based on empirical Bayesian Approach
//' @name getIntensity
//' @param S  Sample covariance matrix.
//' @param target Target matrix.
//' @param lambda A sequence of lambda for shrinkage
//' @param n Number of the samples
//'
//' @return Shrinkage intensity
//' @export

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double getIntensity(arma::mat S, arma::mat target, arma::vec lambda, int n) {
  const int p = S.n_cols;
  const int a = lambda.n_elem;
  arma::vec coef = lambda / (1-lambda);
  arma::vec beta = coef * n + p + 1;
  arma::vec termc(a);
  arma::vec termd(a);
  double val1;
  double sign1;
  double val2;
  double sign2;
  termc.zeros();
  termd.zeros();
  for (int i = 0; i < a; i++){
    for(int j = 0; j < p; j++){
      double c = (beta(i) + (1- j));
      termc(i) += lgamma((c + n) / 2) - lgamma(c / 2);
    }
    log_det(val1, sign1, coef(i) * target);
    termd(i) = beta(i) * val1;
    log_det(val2, sign2, coef(i) * target + S);
    termd(i) -= (beta(i) + n) * val2;
  }
  termd/=2;
  const int index = index_max(-(((n * p) / 2) * log(n * arma::datum::pi)) + termc + termd);
  return lambda[index];
}
