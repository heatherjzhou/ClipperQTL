#include <RcppArmadillo.h>
using namespace arma;

//Given Y, the n*p response matrix,
//and X, the n*K predictor matrix (without a column of ones),
//regress each column of Y against XCombined (with a column of ones),
//and return YResid.t(), the p*n residual matrix.
//[[Rcpp::export()]]
mat residualizeCpp(const mat Y,const mat X){

  // auto start=std::chrono::high_resolution_clock::now();

  int n=X.n_rows; //Sample size.
  vec vecOfOnes(n,fill::ones); //n*1.
  mat XCombined=join_rows(vecOfOnes,X); //n*(K+1).

  mat coef=(XCombined.t()*XCombined).i()*XCombined.t()*Y; //(K+1)*p.
  mat YResid=Y-XCombined*coef; //n*p.

  // auto stop=std::chrono::high_resolution_clock::now();
  // auto duration=std::chrono::duration_cast<std::chrono::milliseconds>(stop-start);
  // cout<<"The C++ code took "<<duration.count()<<" milliseconds.\n";

  return YResid.t(); //p*n.
}



