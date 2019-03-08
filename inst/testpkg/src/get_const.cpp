#include "get_const.h"

extern "C" {

// Get all const columns.

SEXP get_const_all_numeric (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto mat=beachmat::simple_numeric_matrix(in);
    return get_const_all<Rcpp::NumericVector, Rcpp::NumericMatrix>(&mat, order);
    END_RCPP
}

SEXP get_const_all_integer (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto mat=beachmat::simple_integer_matrix(in);
    return get_const_all<Rcpp::IntegerVector, Rcpp::IntegerMatrix>(&mat, order);
    END_RCPP
}

SEXP get_const_all_logical (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto mat=beachmat::simple_logical_matrix(in);
    return get_const_all<Rcpp::LogicalVector, Rcpp::LogicalMatrix>(&mat, order);
    END_RCPP
}

SEXP get_const_all_character (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto mat=beachmat::simple_character_matrix(in);
    return get_const_all<Rcpp::CharacterVector, Rcpp::CharacterMatrix>(&mat, order);
    END_RCPP
}

// Get const column slices.

SEXP get_const_slice_numeric (SEXP in, SEXP order, SEXP bounds) {
    BEGIN_RCPP
    auto mat=beachmat::simple_numeric_matrix(in);
    return get_const_slice<Rcpp::NumericVector, Rcpp::NumericMatrix>(&mat, order, bounds);
    END_RCPP
}

SEXP get_const_slice_integer (SEXP in, SEXP order, SEXP bounds) {
    BEGIN_RCPP
    auto mat=beachmat::simple_integer_matrix(in);
    return get_const_slice<Rcpp::IntegerVector, Rcpp::IntegerMatrix>(&mat, order, bounds);
    END_RCPP
}

SEXP get_const_slice_logical (SEXP in, SEXP order, SEXP bounds) {
    BEGIN_RCPP
    auto mat=beachmat::simple_logical_matrix(in);
    return get_const_slice<Rcpp::LogicalVector, Rcpp::LogicalMatrix>(&mat, order, bounds);
    END_RCPP
}

SEXP get_const_slice_character (SEXP in, SEXP order, SEXP bounds) {
    BEGIN_RCPP
    auto mat=beachmat::simple_character_matrix(in);
    return get_const_slice<Rcpp::CharacterVector, Rcpp::CharacterMatrix>(&mat, order, bounds);
    END_RCPP
}

// Get variable const column slices.

SEXP get_const_varslice_numeric (SEXP in, SEXP order, SEXP bounds) {
    BEGIN_RCPP
    auto mat=beachmat::simple_numeric_matrix(in);
    return get_const_varslice<Rcpp::NumericVector>(&mat, order, bounds);
    END_RCPP
}

SEXP get_const_varslice_integer (SEXP in, SEXP order, SEXP bounds) {
    BEGIN_RCPP
    auto mat=beachmat::simple_integer_matrix(in);
    return get_const_varslice<Rcpp::IntegerVector>(&mat, order, bounds);
    END_RCPP
}

SEXP get_const_varslice_logical (SEXP in, SEXP order, SEXP bounds) {
    BEGIN_RCPP
    auto mat=beachmat::simple_logical_matrix(in);
    return get_const_varslice<Rcpp::LogicalVector>(&mat, order, bounds);
    END_RCPP
}

SEXP get_const_varslice_character (SEXP in, SEXP order, SEXP bounds) {
    BEGIN_RCPP
    auto mat=beachmat::simple_character_matrix(in);
    return get_const_varslice<Rcpp::CharacterVector>(&mat, order, bounds);
    END_RCPP
}

}
