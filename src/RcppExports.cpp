// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// getDataGenotypeCpp
Rcpp::List getDataGenotypeCpp(const std::string chrCurr, const arma::vec geneTSSs, const std::string genotypeFile, const std::string tabixProgram, const std::string tempFileName, const arma::vec sampleIndices, const int cisDistance, const double MAFThreshold, const int MASamplesThreshold);
RcppExport SEXP _ClipperQTL_getDataGenotypeCpp(SEXP chrCurrSEXP, SEXP geneTSSsSEXP, SEXP genotypeFileSEXP, SEXP tabixProgramSEXP, SEXP tempFileNameSEXP, SEXP sampleIndicesSEXP, SEXP cisDistanceSEXP, SEXP MAFThresholdSEXP, SEXP MASamplesThresholdSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const std::string >::type chrCurr(chrCurrSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type geneTSSs(geneTSSsSEXP);
    Rcpp::traits::input_parameter< const std::string >::type genotypeFile(genotypeFileSEXP);
    Rcpp::traits::input_parameter< const std::string >::type tabixProgram(tabixProgramSEXP);
    Rcpp::traits::input_parameter< const std::string >::type tempFileName(tempFileNameSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type sampleIndices(sampleIndicesSEXP);
    Rcpp::traits::input_parameter< const int >::type cisDistance(cisDistanceSEXP);
    Rcpp::traits::input_parameter< const double >::type MAFThreshold(MAFThresholdSEXP);
    Rcpp::traits::input_parameter< const int >::type MASamplesThreshold(MASamplesThresholdSEXP);
    rcpp_result_gen = Rcpp::wrap(getDataGenotypeCpp(chrCurr, geneTSSs, genotypeFile, tabixProgram, tempFileName, sampleIndices, cisDistance, MAFThreshold, MASamplesThreshold));
    return rcpp_result_gen;
END_RCPP
}
// getTableMaxAbsCorsCpp
arma::mat getTableMaxAbsCorsCpp(const arma::mat Y, const arma::mat XRaw, const arma::mat dataCovariates, const arma::vec geneTSSs, const arma::vec SNPPositions, const std::string approach, const int B);
RcppExport SEXP _ClipperQTL_getTableMaxAbsCorsCpp(SEXP YSEXP, SEXP XRawSEXP, SEXP dataCovariatesSEXP, SEXP geneTSSsSEXP, SEXP SNPPositionsSEXP, SEXP approachSEXP, SEXP BSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type XRaw(XRawSEXP);
    Rcpp::traits::input_parameter< const arma::mat >::type dataCovariates(dataCovariatesSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type geneTSSs(geneTSSsSEXP);
    Rcpp::traits::input_parameter< const arma::vec >::type SNPPositions(SNPPositionsSEXP);
    Rcpp::traits::input_parameter< const std::string >::type approach(approachSEXP);
    Rcpp::traits::input_parameter< const int >::type B(BSEXP);
    rcpp_result_gen = Rcpp::wrap(getTableMaxAbsCorsCpp(Y, XRaw, dataCovariates, geneTSSs, SNPPositions, approach, B));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_ClipperQTL_getDataGenotypeCpp", (DL_FUNC) &_ClipperQTL_getDataGenotypeCpp, 9},
    {"_ClipperQTL_getTableMaxAbsCorsCpp", (DL_FUNC) &_ClipperQTL_getTableMaxAbsCorsCpp, 7},
    {NULL, NULL, 0}
};

RcppExport void R_init_ClipperQTL(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
