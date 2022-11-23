// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// myFastLM
arma::vec myFastLM(const arma::mat X, const arma::vec y);
RcppExport SEXP _ClipperQTL_myFastLM(SEXP XSEXP, SEXP ySEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type y(ySEXP);
    rcpp_result_gen = Rcpp::wrap(myFastLM(X, y));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ClipperQTL_myFastLM", (DL_FUNC) &_ClipperQTL_myFastLM, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_ClipperQTL(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
