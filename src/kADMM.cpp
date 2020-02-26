//' ADMM algorithm for estimating sparse precision matrix based on prior knowledge
//' @name kADMM
//' @param Sigma Covariance matrix.
//' @param gammav Regularization term.
//' @param rho Penalty parameter of the augmented Lagrangian function.
//' @param max_iter Maximum number of the iterations.
//' @param numvariable Number of the variables (e.g., Transcription factors)
//' @param prior A zero-one matrix that denotes available prior knowledge
//'
//' @return A Sparse precision matrix
//' @export

#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace arma;
// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
arma::mat kADMM(arma::mat& Sigma, const double gammav, const double rho, int max_iter, const double tol, int numvariable, arma::mat& prior){

  const int abstol = static_cast<double>(2) * tol;
  const int p = Sigma.n_cols;
  double a   = 0.0;
  double b   = 0.0;
  double eps_primal   = 0.0;
  double eps_dual   = 0.0;


  arma::mat Theta(p,p,fill::zeros);
  arma::mat Z(p,p,fill::zeros);
  arma::mat Zold(p,p,fill::zeros);
  arma::mat U(p,p,fill::zeros);

  arma::colvec es(p,fill::zeros);
  arma::colvec thetai(p,fill::zeros);
  arma::mat Q(p,p,fill::zeros);
  double pgammav;

  for (int k = 0;k < max_iter;k++){
    cout << "KDGN : (ADMM) Iteration " << k+1 << " of " << max_iter << " is running. \n";

    //Theta step
    eig_sym(es,Q,rho*(Z-U)-Sigma);
    for (int i = 0;i < p;i++){
      thetai(i) = (es(i)+sqrt(pow(es(i),2)+4*rho))/(2*rho);
    }
    Theta = Q*arma::diagmat(thetai)*Q.t();

    //Z step
    for (int i=0;i<numvariable;i++){
      for (int j=0;j<p;j++){
        double soft1 = 0;
        double soft2 = 0;
        if(prior(i,j)!=0) {pgammav = 0;}
        else{pgammav = gammav/rho; }
        if (Theta(i,j)+U(i,j) > pgammav){soft1 = Theta(i,j) + U(i,j) - pgammav;}
        if (Theta(i,j)+U(i,j) < -pgammav){soft2 = -(Theta(i,j) + U(i,j)) - pgammav;}
        Z(i,j) = soft1 - soft2;
        Z(j,i) = Z(i,j);
      }
    }

    //U step
    U = U + (Theta - Z);

    //Convergence Test
    a = arma::norm(Theta - Z, "fro");
    b = arma::norm(-rho*(Z-Zold),"fro");
    if (norm(Theta,"fro")>norm(Z,"fro")){
      eps_primal = static_cast<double>(p)*abstol + tol*norm(Theta,"fro");
    }else {
      eps_primal = static_cast<double>(p)*abstol + tol*norm(Z,"fro");
    }
    eps_dual = static_cast<double>(p)*abstol + tol*norm(rho*U, "fro");
    if ((a<eps_primal)&&(b<eps_dual)){
      break;
    }
  }
  cout << "KDGN: Precision matrix is returned.\n";
  return (Theta);
}
