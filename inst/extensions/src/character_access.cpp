#include "character_access.h"
#include "aaron_matrix.h"

// Demonstrating with integer matrices.

typedef AaronMatrix<Rcpp::String, Rcpp::StringVector, Rcpp::StringMatrix> AaronStrMat;

// Constructor, destructors and clones.

void * create_character (SEXP incoming) {
    return static_cast<void*>(new AaronStrMat(incoming));
}

void destroy_character(void * ptr) {
    delete static_cast<AaronStrMat*>(ptr);
    return;
}

void * clone_character(void * ptr) {
    AaronStrMat* old=static_cast<AaronStrMat*>(ptr);
    return static_cast<void*>(new AaronStrMat(*old));
}

// Basic getters

void get_dim_character(void* ptr, size_t& nr, size_t& nc){ 
    AaronStrMat* thing=static_cast<AaronStrMat*>(ptr);
    nr=thing->get_nrow();
    nc=thing->get_ncol();
    return;
}

Rcpp::String load_character(void * ptr, size_t r, size_t c) {
    return static_cast<AaronStrMat*>(ptr)->get(r, c);	
}

void load_row_character(void * ptr, size_t r, Rcpp::StringVector::iterator out, size_t first, size_t last) {
    return static_cast<AaronStrMat*>(ptr)->get_row(r, out, first, last);
}

void load_col_character(void * ptr, size_t c, Rcpp::StringVector::iterator out, size_t first, size_t last) {
    return static_cast<AaronStrMat*>(ptr)->get_col(c, out, first, last);
}

// Special getters

Rcpp::StringVector::iterator load_const_col_character(void* ptr, size_t c, Rcpp::StringVector::iterator it, size_t first, size_t last) {
    return static_cast<AaronStrMat*>(ptr)->get_const_col(c, first, last);
}

size_t load_const_col_indexed_character(void* ptr, size_t c, Rcpp::IntegerVector::iterator& iIt, Rcpp::StringVector::iterator& vIt, size_t first, size_t last) {
    return static_cast<AaronStrMat*>(ptr)->get_const_col_indexed(c, iIt, vIt, first, last);
}

// Multi getters

void load_rows_character(void * ptr, Rcpp::IntegerVector::iterator r, size_t n, Rcpp::StringVector::iterator out, size_t first, size_t last) {
    return static_cast<AaronStrMat*>(ptr)->get_rows(r, n, out, first, last);
}

void load_cols_character(void * ptr, Rcpp::IntegerVector::iterator c, size_t n, Rcpp::StringVector::iterator out, size_t first, size_t last) {
    return static_cast<AaronStrMat*>(ptr)->get_cols(c, n, out, first, last);
}

