//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
using namespace arma;

//[[Rcpp::export()]]
vec testFunction1(const mat X){

  vec toReturn=max(X,1); //1 means to return statistic for each row.

  return toReturn;
}

//[[Rcpp::export()]]
double testFunction2(const vec x){

  double toReturn=max(x);

  return toReturn;
}





