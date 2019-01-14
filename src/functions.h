#ifndef BEACHMAT_FUNCTIONS_H
#define BEACHMAT_FUNCTIONS_H

#include "Rcpp.h"

extern "C" { 

SEXP rechunk_matrix(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);

SEXP find_chunks(SEXP);

}

#endif
