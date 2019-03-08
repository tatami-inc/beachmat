#include "Rcpp.h"

extern "C" {

void * AaronMatrix_character_input_create(SEXP);

void AaronMatrix_character_input_destroy(void *);

void * AaronMatrix_character_input_clone(void *);

void AaronMatrix_character_input_dim(void*, size_t*, size_t*);
 
void AaronMatrix_character_input_get(void *, size_t, size_t, Rcpp::String*);

void AaronMatrix_character_input_getRow(void *, size_t, Rcpp::StringVector::iterator*, size_t, size_t);

void AaronMatrix_character_input_getCol(void *, size_t, Rcpp::StringVector::iterator*, size_t, size_t);

void AaronMatrix_character_input_getConstCol(void*, size_t, Rcpp::StringVector::iterator*, size_t, size_t, Rcpp::StringVector::iterator*);

size_t AaronMatrix_character_input_getConstColIndexed(void*, size_t, Rcpp::IntegerVector::iterator*, Rcpp::StringVector::iterator*, size_t, size_t);

void AaronMatrix_character_input_getRows(void *, Rcpp::IntegerVector::iterator*, size_t, Rcpp::StringVector::iterator*, size_t, size_t);

void AaronMatrix_character_input_getCols(void *, Rcpp::IntegerVector::iterator*, size_t, Rcpp::StringVector::iterator*, size_t, size_t);

void * AaronMatrix_integer_input_create (SEXP);

void AaronMatrix_integer_input_destroy (void *);

void * AaronMatrix_integer_input_clone (void *);

void AaronMatrix_integer_input_dim(void*, size_t*, size_t*);
 
void AaronMatrix_integer_input_get(void *, size_t, size_t, int*);

void AaronMatrix_integer_input_getRow_integer(void *, size_t, Rcpp::IntegerVector::iterator*, size_t, size_t);

void AaronMatrix_integer_input_getCol_integer(void *, size_t, Rcpp::IntegerVector::iterator*, size_t, size_t);

void AaronMatrix_integer_input_getRow_numeric(void *, size_t, Rcpp::NumericVector::iterator*, size_t, size_t);

void AaronMatrix_integer_input_getCol_numeric(void *, size_t, Rcpp::NumericVector::iterator*, size_t, size_t);

void AaronMatrix_integer_input_getConstCol(void*, size_t, Rcpp::IntegerVector::iterator*, size_t, size_t, Rcpp::IntegerVector::iterator*);

size_t AaronMatrix_integer_input_getConstColIndexed(void*, size_t, Rcpp::IntegerVector::iterator*, Rcpp::IntegerVector::iterator*, size_t, size_t);

void AaronMatrix_integer_input_getRows_integer(void *, Rcpp::IntegerVector::iterator*, size_t, Rcpp::IntegerVector::iterator*, size_t, size_t);

void AaronMatrix_integer_input_getCols_integer(void *, Rcpp::IntegerVector::iterator*, size_t, Rcpp::IntegerVector::iterator*, size_t, size_t);

void AaronMatrix_integer_input_getRows_numeric(void *, Rcpp::IntegerVector::iterator*, size_t, Rcpp::NumericVector::iterator*, size_t, size_t);

void AaronMatrix_integer_input_getCols_numeric(void *, Rcpp::IntegerVector::iterator*, size_t, Rcpp::NumericVector::iterator*, size_t, size_t);

}
