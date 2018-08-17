#include "Rcpp.h"

extern "C" {

void * create_integer (SEXP);

void destroy_integer(void *);

void * clone_integer(void *);

// Basic getters.

void get_dim_integer(void*, size_t&, size_t&);

int load_integer(void*, size_t, size_t);

void load_row2int_integer(void *, size_t, Rcpp::IntegerVector::iterator, size_t, size_t);

void load_col2int_integer(void *, size_t, Rcpp::IntegerVector::iterator, size_t, size_t);

void load_row2dbl_integer(void *, size_t, Rcpp::NumericVector::iterator, size_t, size_t);

void load_col2dbl_integer(void *, size_t, Rcpp::NumericVector::iterator, size_t, size_t);

// Special getters.

Rcpp::IntegerVector::iterator load_const_col_integer(void *, size_t, Rcpp::IntegerVector::iterator, size_t, size_t);

size_t load_const_col_indexed_integer(void *, size_t, Rcpp::IntegerVector::iterator&, Rcpp::IntegerVector::iterator&, size_t, size_t);

// Multi getters.

void load_rows2int_integer(void *, Rcpp::IntegerVector::iterator, size_t, Rcpp::IntegerVector::iterator, size_t, size_t);

void load_cols2int_integer(void *, Rcpp::IntegerVector::iterator, size_t, Rcpp::IntegerVector::iterator, size_t, size_t);

void load_rows2dbl_integer(void *, Rcpp::IntegerVector::iterator, size_t, Rcpp::NumericVector::iterator, size_t, size_t);

void load_cols2dbl_integer(void *, Rcpp::IntegerVector::iterator, size_t, Rcpp::NumericVector::iterator, size_t, size_t);

}
