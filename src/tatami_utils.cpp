#include "tatamize.h"
#include "Rcpp.h"

//[[Rcpp::export(rng=false)]]
Rcpp::IntegerVector tatami_dim(SEXP input) {
    auto shared = extract_NumericMatrix_shared(input);
    return Rcpp::IntegerVector::create(shared->nrow(), shared->ncol());
}

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector tatami_column(SEXP input, int i) {
    auto shared = extract_NumericMatrix_shared(input);
    Rcpp::NumericVector output(shared->nrow());
    auto wrk = shared->dense_column();
    wrk->fetch_copy(i-1, static_cast<double*>(output.begin()));
    return output;
}

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector tatami_row(SEXP input, int i) {
    auto shared = extract_NumericMatrix_shared(input);
    Rcpp::NumericVector output(shared->ncol());
    auto wrk = shared->dense_row();
    wrk->fetch_copy(i-1, static_cast<double*>(output.begin()));
    return output;
}
