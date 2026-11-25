#include "Rtatami.h"
#include "Rcpp.h"

#include "sanisizer/sanisizer.hpp"
#include "tatami/tatami.hpp"

#include <vector>
#include <cstddef>
#include <stdexcept>

#include "utils.h"

//[[Rcpp::export(rng=false)]]
Rcpp::IntegerVector tatami_dim(SEXP raw_input) {
    Rtatami::BoundNumericPointer input(raw_input);
    const auto& shared = input->ptr;
    return Rcpp::IntegerVector::create(shared->nrow(), shared->ncol());
}

//[[Rcpp::export(rng=false)]]
Rcpp::LogicalVector tatami_is_sparse(SEXP raw_input) {
    Rtatami::BoundNumericPointer input(raw_input);
    const auto& shared = input->ptr;
    return Rcpp::LogicalVector::create(shared->is_sparse());
}

//[[Rcpp::export(rng=false)]]
Rcpp::LogicalVector tatami_prefer_rows(SEXP raw_input) {
    Rtatami::BoundNumericPointer input(raw_input);
    const auto& shared = input->ptr;
    return Rcpp::LogicalVector::create(shared->prefer_rows());
}

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector tatami_get(SEXP raw_input, int i, bool row) {
    Rtatami::BoundNumericPointer input(raw_input);
    const auto NR = input->ptr->nrow();
    const auto NC = input->ptr->ncol();

    const auto limit = (row ? NR : NC);
    if (i < 1 || i > limit) {
        throw std::runtime_error("'i' is out of range");
    }

    auto output = sanisizer::create<Rcpp::NumericVector>(row ? NC : NR);
    auto optr = static_cast<double*>(output.begin());

    auto wrk = tatami::consecutive_extractor<false>(*(input->ptr), row, i - 1, 1);
    auto ptr = wrk->fetch(optr);
    tatami::copy_n(ptr, output.size(), optr);
    return output;
}
