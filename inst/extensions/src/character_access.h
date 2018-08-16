#include "Rcpp.h"

extern "C" {

void * create_character (SEXP incoming);

void destroy_character(void *);

void * clone_character(void *);

// Basic getters

void get_dim_character(void*, size_t*, size_t*);

Rcpp::String load_character(void *, size_t, size_t);

void load_row_character(void *, size_t, Rcpp::StringVector::iterator, size_t, size_t);

void load_col_character(void *, size_t, Rcpp::StringVector::iterator, size_t, size_t);

// Multi getters

void load_rows_character(void *, Rcpp::IntegerVector::iterator, size_t, Rcpp::StringVector::iterator, size_t, size_t);

void load_cols_character(void *, Rcpp::IntegerVector::iterator, size_t, Rcpp::StringVector::iterator, size_t, size_t);

}
