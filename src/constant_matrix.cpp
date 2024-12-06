#include "Rtatami.h"

//[[Rcpp::export(rng=false)]]
SEXP initialize_constant_matrix(int nrow, int ncol, double val) {
    auto output = Rtatami::new_BoundNumericMatrix();
    output->ptr.reset(new tatami::ConstantMatrix(nrow, ncol, val));
    return output;
}
