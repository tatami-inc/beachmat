#include "integer_access.h"
#include "aaron_matrix.h"

// Demonstrating with integer matrices.

typedef AaronMatrix<int, Rcpp::IntegerVector, Rcpp::IntegerMatrix> AaronIntMat;

// Constructor, destructors and clones.

void * create_integer (SEXP incoming) {
    return static_cast<void*>(new AaronIntMat(incoming));
}

void destroy_integer(void * ptr) {
    delete static_cast<AaronIntMat*>(ptr);
    return;
}

void * clone_integer(void * ptr) {
    AaronIntMat* old=static_cast<AaronIntMat*>(ptr);
    return static_cast<void*>(new AaronIntMat(*old));
}

// Basic getters

void get_dim_integer(void* ptr, size_t* nr, size_t* nc){ 
    AaronIntMat* thing=static_cast<AaronIntMat*>(ptr);
    *nr=thing->get_nrow();
    *nc=thing->get_ncol();
    return;
}

int load_integer(void * ptr, size_t r, size_t c) {
    return static_cast<AaronIntMat*>(ptr)->get(r, c);	
}

void load_row2int_integer(void * ptr, size_t r, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
    return static_cast<AaronIntMat*>(ptr)->get_row(r, out, first, last);
}

void load_col2int_integer(void * ptr, size_t c, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
    return static_cast<AaronIntMat*>(ptr)->get_col(c, out, first, last);
}

void load_row2dbl_integer(void * ptr, size_t r, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
    return static_cast<AaronIntMat*>(ptr)->get_row(r, out, first, last);
}

void load_col2dbl_integer(void * ptr, size_t c, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
    return static_cast<AaronIntMat*>(ptr)->get_col(c, out, first, last);
}

// Multi getters

void load_rows2int_integer(void * ptr, Rcpp::IntegerVector::iterator r, size_t n, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
    return static_cast<AaronIntMat*>(ptr)->get_rows(r, n, out, first, last);
}

void load_cols2int_integer(void * ptr, Rcpp::IntegerVector::iterator c, size_t n, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
    return static_cast<AaronIntMat*>(ptr)->get_cols(c, n, out, first, last);
}

void load_rows2dbl_integer(void * ptr, Rcpp::IntegerVector::iterator r, size_t n, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
    return static_cast<AaronIntMat*>(ptr)->get_rows(r, n, out, first, last);
}

void load_cols2dbl_integer(void * ptr, Rcpp::IntegerVector::iterator c, size_t n, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
    return static_cast<AaronIntMat*>(ptr)->get_cols(c, n, out, first, last);
}
