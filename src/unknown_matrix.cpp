#include "Rtatami.h"
#include "Rcpp.h"
#include "tatami_r/tatami_r.hpp"

//[[Rcpp::export(rng=false)]]
SEXP initialize_unknown_matrix(Rcpp::RObject input) {
    auto output = Rtatami::new_BoundNumericMatrix();
    output->original = input;
    output->ptr.reset(new tatami_r::UnknownMatrix<double, int>(input));
    return output;
}
