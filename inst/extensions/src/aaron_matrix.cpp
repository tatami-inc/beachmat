#include "aaron_matrix.h"

// Demonstrating with integer matrices only, for simplicity.

typedef AaronMatrix<int, Rcpp::IntegerVector, Rcpp::IntegerMatrix> AaronIntMat;

extern "C" {

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

void get_dim_integer(void* ptr, int* nr, int* nc){ 
    AaronIntMat* thing=static_cast<AaronIntMat*>(ptr);
    *nr=thing->get_nrow();
    *nc=thing->get_ncol();
    return;
}

int load_integer(void * ptr, int r, int c) {
    return static_cast<AaronIntMat*>(ptr)->get(r, c);	
}

void load_row2int_integer(void * ptr, int r, int* out, int first, int last) {
    return static_cast<AaronIntMat*>(ptr)->get_row(r, out, first, last);
}

void load_col2int_integer(void * ptr, int c, int* out, int first, int last) {
    return static_cast<AaronIntMat*>(ptr)->get_col(c, out, first, last);
}

void load_row2dbl_integer(void * ptr, int r, double* out, int first, int last) {
    return static_cast<AaronIntMat*>(ptr)->get_row(r, out, first, last);
}

void load_col2dbl_integer(void * ptr, int c, double* out, int first, int last) {
    return static_cast<AaronIntMat*>(ptr)->get_col(c, out, first, last);
}

// Multi getters

void load_rows2int_integer(void * ptr, int* r, int n, int* out, int first, int last) {
    return static_cast<AaronIntMat*>(ptr)->get_rows(r, n, out, first, last);
}

void load_cols2int_integer(void * ptr, int* c, int n, int* out, int first, int last) {
    return static_cast<AaronIntMat*>(ptr)->get_cols(c, n, out, first, last);
}

void load_rows2dbl_integer(void * ptr, int* r, int n, double* out, int first, int last) {
    return static_cast<AaronIntMat*>(ptr)->get_rows(r, n, out, first, last);
}

void load_cols2dbl_integer(void * ptr, int* c, int n, double* out, int first, int last) {
    return static_cast<AaronIntMat*>(ptr)->get_cols(c, n, out, first, last);
}

}

#include "R_ext/Rdynload.h"

#define REGISTER(x) R_RegisterCCallable("morebeach", #x, reinterpret_cast<DL_FUNC>(x))

extern "C" {

// Defining the initialization function.

void R_init_morebeach(DllInfo *info) {
    REGISTER(create_integer);
    REGISTER(destroy_integer);
    REGISTER(clone_integer);

    REGISTER(get_dim_integer);
    REGISTER(load_integer);
    REGISTER(load_row2int_integer);
    REGISTER(load_col2int_integer);
    REGISTER(load_row2dbl_integer);
    REGISTER(load_col2dbl_integer);

    REGISTER(load_rows2int_integer);
    REGISTER(load_cols2int_integer);
    REGISTER(load_rows2dbl_integer);
    REGISTER(load_cols2dbl_integer);
}

}
