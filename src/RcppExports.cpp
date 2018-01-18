// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <Rcpp.h>

using namespace Rcpp;

// distance
Rcpp::NumericMatrix distance(Rcpp::NumericMatrix start, Rcpp::NumericMatrix end);
RcppExport SEXP _vegdistmod_distance(SEXP startSEXP, SEXP endSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type start(startSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type end(endSEXP);
    rcpp_result_gen = Rcpp::wrap(distance(start, end));
    return rcpp_result_gen;
END_RCPP
}
// findcoord
Rcpp::NumericVector findcoord(double lon, double lat, double dist, double brng);
RcppExport SEXP _vegdistmod_findcoord(SEXP lonSEXP, SEXP latSEXP, SEXP distSEXP, SEXP brngSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< double >::type lon(lonSEXP);
    Rcpp::traits::input_parameter< double >::type lat(latSEXP);
    Rcpp::traits::input_parameter< double >::type dist(distSEXP);
    Rcpp::traits::input_parameter< double >::type brng(brngSEXP);
    rcpp_result_gen = Rcpp::wrap(findcoord(lon, lat, dist, brng));
    return rcpp_result_gen;
END_RCPP
}
// latlonfromcell
Rcpp::NumericMatrix latlonfromcell(Rcpp::NumericVector cells, Rcpp::NumericVector extent, int nrow, int ncol);
RcppExport SEXP _vegdistmod_latlonfromcell(SEXP cellsSEXP, SEXP extentSEXP, SEXP nrowSEXP, SEXP ncolSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type cells(cellsSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericVector >::type extent(extentSEXP);
    Rcpp::traits::input_parameter< int >::type nrow(nrowSEXP);
    Rcpp::traits::input_parameter< int >::type ncol(ncolSEXP);
    rcpp_result_gen = Rcpp::wrap(latlonfromcell(cells, extent, nrow, ncol));
    return rcpp_result_gen;
END_RCPP
}
// pythagorean
Rcpp::NumericMatrix pythagorean(Rcpp::NumericMatrix start, Rcpp::NumericMatrix end);
RcppExport SEXP _vegdistmod_pythagorean(SEXP startSEXP, SEXP endSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type start(startSEXP);
    Rcpp::traits::input_parameter< Rcpp::NumericMatrix >::type end(endSEXP);
    rcpp_result_gen = Rcpp::wrap(pythagorean(start, end));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_vegdistmod_distance", (DL_FUNC) &_vegdistmod_distance, 2},
    {"_vegdistmod_findcoord", (DL_FUNC) &_vegdistmod_findcoord, 4},
    {"_vegdistmod_latlonfromcell", (DL_FUNC) &_vegdistmod_latlonfromcell, 4},
    {"_vegdistmod_pythagorean", (DL_FUNC) &_vegdistmod_pythagorean, 2},
    {NULL, NULL, 0}
};

RcppExport void R_init_vegdistmod(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
