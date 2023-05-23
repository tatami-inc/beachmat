#include "Rtatami.h"
#include "Rcpp.h"

//[[Rcpp::export(rng=false)]]
Rcpp::IntegerVector tatami_dim(SEXP raw_input) {
    Rtatami::BoundNumericPointer input(raw_input);
    const auto& shared = input->ptr;
    return Rcpp::IntegerVector::create(shared->nrow(), shared->ncol());
}

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector tatami_column(SEXP raw_input, int i) {
    Rtatami::BoundNumericPointer input(raw_input);
    const auto& shared = input->ptr;
    Rcpp::NumericVector output(shared->nrow());
    auto wrk = shared->dense_column();
    wrk->fetch_copy(i-1, static_cast<double*>(output.begin()));
    return output;
}

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector tatami_row(SEXP raw_input, int i) {
    Rtatami::BoundNumericPointer input(raw_input);
    const auto& shared = input->ptr;
    Rcpp::NumericVector output(shared->ncol());
    auto wrk = shared->dense_row();
    wrk->fetch_copy(i-1, static_cast<double*>(output.begin()));
    return output;
}
