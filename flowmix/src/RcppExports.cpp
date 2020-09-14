// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// soft_threshC
NumericVector soft_threshC(const NumericVector& a, const double& b);
RcppExport SEXP _flowmix_soft_threshC(SEXP aSEXP, SEXP bSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type a(aSEXP);
    Rcpp::traits::input_parameter< const double& >::type b(bSEXP);
    rcpp_result_gen = Rcpp::wrap(soft_threshC(a, b));
    return rcpp_result_gen;
END_RCPP
}
// wvec_updateC
NumericVector wvec_updateC(const NumericVector& b1, const NumericVector& uw, const double& lambda, const double& rho);
RcppExport SEXP _flowmix_wvec_updateC(SEXP b1SEXP, SEXP uwSEXP, SEXP lambdaSEXP, SEXP rhoSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type b1(b1SEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type uw(uwSEXP);
    Rcpp::traits::input_parameter< const double& >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const double& >::type rho(rhoSEXP);
    rcpp_result_gen = Rcpp::wrap(wvec_updateC(b1, uw, lambda, rho));
    return rcpp_result_gen;
END_RCPP
}
// rowSumsC
NumericVector rowSumsC(NumericMatrix x);
RcppExport SEXP _flowmix_rowSumsC(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< NumericMatrix >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rowSumsC(x));
    return rcpp_result_gen;
END_RCPP
}
// rowSumsC_arma
NumericVector rowSumsC_arma(const arma::mat& x);
RcppExport SEXP _flowmix_rowSumsC_arma(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rowSumsC_arma(x));
    return rcpp_result_gen;
END_RCPP
}
// rowSumsC2_arma
NumericVector rowSumsC2_arma(const arma::mat& x);
RcppExport SEXP _flowmix_rowSumsC2_arma(SEXP xSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type x(xSEXP);
    rcpp_result_gen = Rcpp::wrap(rowSumsC2_arma(x));
    return rcpp_result_gen;
END_RCPP
}
// projCmatC
arma::mat projCmatC(const arma::mat& mat, const double& C);
RcppExport SEXP _flowmix_projCmatC(SEXP matSEXP, SEXP CSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< const double& >::type C(CSEXP);
    rcpp_result_gen = Rcpp::wrap(projCmatC(mat, C));
    return rcpp_result_gen;
END_RCPP
}
// Z_updateC
arma::mat Z_updateC(arma::mat Xbeta1, arma::mat Uz, double C, double rho, int dimdat, int TT);
RcppExport SEXP _flowmix_Z_updateC(SEXP Xbeta1SEXP, SEXP UzSEXP, SEXP CSEXP, SEXP rhoSEXP, SEXP dimdatSEXP, SEXP TTSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type Xbeta1(Xbeta1SEXP);
    Rcpp::traits::input_parameter< arma::mat >::type Uz(UzSEXP);
    Rcpp::traits::input_parameter< double >::type C(CSEXP);
    Rcpp::traits::input_parameter< double >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< int >::type dimdat(dimdatSEXP);
    Rcpp::traits::input_parameter< int >::type TT(TTSEXP);
    rcpp_result_gen = Rcpp::wrap(Z_updateC(Xbeta1, Uz, C, rho, dimdat, TT));
    return rcpp_result_gen;
END_RCPP
}
// mv_mult
arma::mat mv_mult(const arma::mat& lhs, const arma::vec& rhs);
RcppExport SEXP _flowmix_mv_mult(SEXP lhsSEXP, SEXP rhsSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type lhs(lhsSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type rhs(rhsSEXP);
    rcpp_result_gen = Rcpp::wrap(mv_mult(lhs, rhs));
    return rcpp_result_gen;
END_RCPP
}
// b_updateC
arma::vec b_updateC(const NumericVector& wvec, const NumericVector& uw, const double& rho, const NumericVector& cvec3_el, const NumericVector& yvec, const arma::mat& DtDinv, const arma::mat& D, const double& N);
RcppExport SEXP _flowmix_b_updateC(SEXP wvecSEXP, SEXP uwSEXP, SEXP rhoSEXP, SEXP cvec3_elSEXP, SEXP yvecSEXP, SEXP DtDinvSEXP, SEXP DSEXP, SEXP NSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type wvec(wvecSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type uw(uwSEXP);
    Rcpp::traits::input_parameter< const double& >::type rho(rhoSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type cvec3_el(cvec3_elSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type yvec(yvecSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type DtDinv(DtDinvSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type D(DSEXP);
    Rcpp::traits::input_parameter< const double& >::type N(NSEXP);
    rcpp_result_gen = Rcpp::wrap(b_updateC(wvec, uw, rho, cvec3_el, yvec, DtDinv, D, N));
    return rcpp_result_gen;
END_RCPP
}
// subtractC2
arma::mat subtractC2(const arma::vec& wt, const arma::mat& mat, const arma::vec& vec, const arma::mat& mat2);
RcppExport SEXP _flowmix_subtractC2(SEXP wtSEXP, SEXP matSEXP, SEXP vecSEXP, SEXP mat2SEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type wt(wtSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type vec(vecSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type mat2(mat2SEXP);
    rcpp_result_gen = Rcpp::wrap(subtractC2(wt, mat, vec, mat2));
    return rcpp_result_gen;
END_RCPP
}
// subtractC3
arma::mat subtractC3(const arma::vec& wt, const arma::mat& mat, const arma::vec& vec);
RcppExport SEXP _flowmix_subtractC3(SEXP wtSEXP, SEXP matSEXP, SEXP vecSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type wt(wtSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type mat(matSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type vec(vecSEXP);
    rcpp_result_gen = Rcpp::wrap(subtractC3(wt, mat, vec));
    return rcpp_result_gen;
END_RCPP
}
// dothisC
arma::mat dothisC(const arma::vec& longwt, const arma::mat& ylong, const arma::mat& mumat, const arma::mat& sigma_half);
RcppExport SEXP _flowmix_dothisC(SEXP longwtSEXP, SEXP ylongSEXP, SEXP mumatSEXP, SEXP sigma_halfSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type longwt(longwtSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type ylong(ylongSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type mumat(mumatSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type sigma_half(sigma_halfSEXP);
    rcpp_result_gen = Rcpp::wrap(dothisC(longwt, ylong, mumat, sigma_half));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_flowmix_soft_threshC", (DL_FUNC) &_flowmix_soft_threshC, 2},
    {"_flowmix_wvec_updateC", (DL_FUNC) &_flowmix_wvec_updateC, 4},
    {"_flowmix_rowSumsC", (DL_FUNC) &_flowmix_rowSumsC, 1},
    {"_flowmix_rowSumsC_arma", (DL_FUNC) &_flowmix_rowSumsC_arma, 1},
    {"_flowmix_rowSumsC2_arma", (DL_FUNC) &_flowmix_rowSumsC2_arma, 1},
    {"_flowmix_projCmatC", (DL_FUNC) &_flowmix_projCmatC, 2},
    {"_flowmix_Z_updateC", (DL_FUNC) &_flowmix_Z_updateC, 6},
    {"_flowmix_mv_mult", (DL_FUNC) &_flowmix_mv_mult, 2},
    {"_flowmix_b_updateC", (DL_FUNC) &_flowmix_b_updateC, 8},
    {"_flowmix_subtractC2", (DL_FUNC) &_flowmix_subtractC2, 4},
    {"_flowmix_subtractC3", (DL_FUNC) &_flowmix_subtractC3, 3},
    {"_flowmix_dothisC", (DL_FUNC) &_flowmix_dothisC, 4},
    {NULL, NULL, 0}
};

RcppExport void R_init_flowmix(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
