#include "tatamize.h"
#include "Rcpp.h"
#include "tatami_r/tatami_r.hpp"

//[[Rcpp::export(rng=false)]]
SEXP initialize_unknown_matrix(Rcpp::RObject input) {
    return new_MatrixChan(new tatami_r::UnknownMatrix<double, int>(input));
}
