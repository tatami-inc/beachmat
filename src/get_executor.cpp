#include "Rtatami.h"

//[[Rcpp::export(rng=false)]]
SEXP get_executor() {
    // Do NOT destruct the global value, hence the false.
    return Rcpp::XPtr<manticore::Executor>(&(tatami_r::executor()), false);
}
