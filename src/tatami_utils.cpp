#include "Rtatami.h"
#include "Rcpp.h"
#include "tatami_stats/tatami_stats.hpp"

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
    auto wrk = shared->dense_column();

    Rcpp::NumericVector output(shared->nrow());
    auto optr = static_cast<double*>(output.begin());
    auto ptr = wrk->fetch(i-1, optr);
    tatami::copy_n(ptr, output.size(), optr);
    return output;
}

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector tatami_row(SEXP raw_input, int i) {
    Rtatami::BoundNumericPointer input(raw_input);
    const auto& shared = input->ptr;
    auto wrk = shared->dense_row();

    Rcpp::NumericVector output(shared->ncol());
    auto optr = static_cast<double*>(output.begin());
    auto ptr = wrk->fetch(i-1, optr);
    tatami::copy_n(ptr, output.size(), optr);
    return output;
}

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector tatami_row_sums(SEXP raw_input, int threads) {
    Rtatami::BoundNumericPointer input(raw_input);
    tatami_stats::sums::Options opt;
    opt.num_threads = threads;
    auto rs = tatami_stats::sums::by_row(input->ptr.get(), opt);
    return Rcpp::NumericVector(rs.begin(), rs.end());
}

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector tatami_column_sums(SEXP raw_input, int threads) {
    Rtatami::BoundNumericPointer input(raw_input);
    tatami_stats::sums::Options opt;
    opt.num_threads = threads;
    auto rs = tatami_stats::sums::by_column(input->ptr.get(), opt);
    return Rcpp::NumericVector(rs.begin(), rs.end());
}
