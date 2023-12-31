// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// CMC_C
List CMC_C(int nDim_input, NumericVector M_input, NumericVector Margin_input, NumericVector Y_input, NumericVector item_used_input);
RcppExport SEXP _CMC_CMC_C(SEXP nDim_inputSEXP, SEXP M_inputSEXP, SEXP Margin_inputSEXP, SEXP Y_inputSEXP, SEXP item_used_inputSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< int >::type nDim_input(nDim_inputSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type M_input(M_inputSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Margin_input(Margin_inputSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type Y_input(Y_inputSEXP);
    Rcpp::traits::input_parameter< NumericVector >::type item_used_input(item_used_inputSEXP);
    rcpp_result_gen = Rcpp::wrap(CMC_C(nDim_input, M_input, Margin_input, Y_input, item_used_input));
    return rcpp_result_gen;
END_RCPP
}
// rcpp_hello_world
List rcpp_hello_world();
RcppExport SEXP _CMC_rcpp_hello_world() {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    rcpp_result_gen = Rcpp::wrap(rcpp_hello_world());
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_CMC_CMC_C", (DL_FUNC) &_CMC_CMC_C, 5},
    {"_CMC_rcpp_hello_world", (DL_FUNC) &_CMC_rcpp_hello_world, 0},
    {NULL, NULL, 0}
};

RcppExport void R_init_CMC(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
