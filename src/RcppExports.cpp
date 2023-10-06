// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// RMST
SEXP RMST(const arma::colvec status, const arma::colvec time, const bool extend, Rcpp::Nullable<double> tau);
RcppExport SEXP _MRMST_RMST(SEXP statusSEXP, SEXP timeSEXP, SEXP extendSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec >::type status(statusSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type time(timeSEXP);
    Rcpp::traits::input_parameter< const bool >::type extend(extendSEXP);
    Rcpp::traits::input_parameter< Rcpp::Nullable<double> >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(RMST(status, time, extend, tau));
    return rcpp_result_gen;
END_RCPP
}
// CalcPsiRMST
SEXP CalcPsiRMST(const arma::colvec status, const arma::colvec time, const double tau);
RcppExport SEXP _MRMST_CalcPsiRMST(SEXP statusSEXP, SEXP timeSEXP, SEXP tauSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec >::type status(statusSEXP);
    Rcpp::traits::input_parameter< const arma::colvec >::type time(timeSEXP);
    Rcpp::traits::input_parameter< const double >::type tau(tauSEXP);
    rcpp_result_gen = Rcpp::wrap(CalcPsiRMST(status, time, tau));
    return rcpp_result_gen;
END_RCPP
}
// GenPerturb
SEXP GenPerturb(const arma::colvec psi, const int n_boot);
RcppExport SEXP _MRMST_GenPerturb(SEXP psiSEXP, SEXP n_bootSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::colvec >::type psi(psiSEXP);
    Rcpp::traits::input_parameter< const int >::type n_boot(n_bootSEXP);
    rcpp_result_gen = Rcpp::wrap(GenPerturb(psi, n_boot));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_MRMST_RMST", (DL_FUNC) &_MRMST_RMST, 4},
    {"_MRMST_CalcPsiRMST", (DL_FUNC) &_MRMST_CalcPsiRMST, 3},
    {"_MRMST_GenPerturb", (DL_FUNC) &_MRMST_GenPerturb, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_MRMST(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
