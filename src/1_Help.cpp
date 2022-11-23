#include <RcppArmadillo.h>
using namespace arma;

//[[Rcpp::export]]
arma::vec myFastLM(const arma::mat X,const arma::vec y){
  vec coef=solve(X,y);
  return coef;
}
