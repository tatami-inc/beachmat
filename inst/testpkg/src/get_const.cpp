#include "get_const.h"

extern "C" {

// Get all const columns.

SEXP get_const_all_numeric (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(in);
    return get_const_all<Rcpp::NumericVector, Rcpp::NumericMatrix>(ptr.get(), order);
    END_RCPP
}

SEXP get_const_all_integer (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_integer_matrix(in);
    return get_const_all<Rcpp::IntegerVector, Rcpp::IntegerMatrix>(ptr.get(), order);
    END_RCPP
}

SEXP get_const_all_logical (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_logical_matrix(in);
    return get_const_all<Rcpp::LogicalVector, Rcpp::LogicalMatrix>(ptr.get(), order);
    END_RCPP
}

SEXP get_const_all_character (SEXP in, SEXP order) {
    BEGIN_RCPP
    auto ptr=beachmat::create_character_matrix(in);
    return get_const_all<Rcpp::CharacterVector, Rcpp::CharacterMatrix>(ptr.get(), order);
    END_RCPP
}

// Get const column slices.

SEXP get_const_slice_numeric (SEXP in, SEXP order, SEXP bounds) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(in);
    return get_const_slice<Rcpp::NumericVector, Rcpp::NumericMatrix>(ptr.get(), order, bounds);
    END_RCPP
}

SEXP get_const_slice_integer (SEXP in, SEXP order, SEXP bounds) {
    BEGIN_RCPP
    auto ptr=beachmat::create_integer_matrix(in);
    return get_const_slice<Rcpp::IntegerVector, Rcpp::IntegerMatrix>(ptr.get(), order, bounds);
    END_RCPP
}

SEXP get_const_slice_logical (SEXP in, SEXP order, SEXP bounds) {
    BEGIN_RCPP
    auto ptr=beachmat::create_logical_matrix(in);
    return get_const_slice<Rcpp::LogicalVector, Rcpp::LogicalMatrix>(ptr.get(), order, bounds);
    END_RCPP
}

SEXP get_const_slice_character (SEXP in, SEXP order, SEXP bounds) {
    BEGIN_RCPP
    auto ptr=beachmat::create_character_matrix(in);
    return get_const_slice<Rcpp::CharacterVector, Rcpp::CharacterMatrix>(ptr.get(), order, bounds);
    END_RCPP
}

// Get variable const column slices.

SEXP get_const_varslice_numeric (SEXP in, SEXP order, SEXP bounds) {
    BEGIN_RCPP
    auto ptr=beachmat::create_numeric_matrix(in);
    return get_const_varslice<Rcpp::NumericVector>(ptr.get(), order, bounds);
    END_RCPP
}

SEXP get_const_varslice_integer (SEXP in, SEXP order, SEXP bounds) {
    BEGIN_RCPP
    auto ptr=beachmat::create_integer_matrix(in);
    return get_const_varslice<Rcpp::IntegerVector>(ptr.get(), order, bounds);
    END_RCPP
}

SEXP get_const_varslice_logical (SEXP in, SEXP order, SEXP bounds) {
    BEGIN_RCPP
    auto ptr=beachmat::create_logical_matrix(in);
    return get_const_varslice<Rcpp::LogicalVector>(ptr.get(), order, bounds);
    END_RCPP
}

SEXP get_const_varslice_character (SEXP in, SEXP order, SEXP bounds) {
    BEGIN_RCPP
    auto ptr=beachmat::create_character_matrix(in);
    return get_const_varslice<Rcpp::CharacterVector>(ptr.get(), order, bounds);
    END_RCPP
}

}
