#include "Rcpp.h"

extern "C" {

void * create_integer (SEXP);

void destroy_integer(void *);

void * clone_integer(void *);

void get_dim_integer(void*, size_t*, size_t*);

int load_integer(void*, size_t, size_t);

void load_row2int_integer(void *, size_t, Rcpp::IntegerVector::iterator, size_t, size_t);

void load_col2int_integer(void *, size_t, Rcpp::IntegerVector::iterator, size_t, size_t);

void load_row2dbl_integer(void *, size_t, Rcpp::NumericVector::iterator, size_t, size_t);

void load_col2dbl_integer(void *, size_t, Rcpp::NumericVector::iterator, size_t, size_t);

void load_rows2int_integer(void *, Rcpp::IntegerVector::iterator, size_t, Rcpp::IntegerVector::iterator, size_t, size_t);

void load_cols2int_integer(void *, Rcpp::IntegerVector::iterator, size_t, Rcpp::IntegerVector::iterator, size_t, size_t);

void load_rows2dbl_integer(void *, Rcpp::IntegerVector::iterator, size_t, Rcpp::NumericVector::iterator, size_t, size_t);

void load_cols2dbl_integer(void *, Rcpp::IntegerVector::iterator, size_t, Rcpp::NumericVector::iterator, size_t, size_t);

}
