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

//[[Rcpp::export(rng=false)]]
SEXP tatami_realize(SEXP raw_input, int threads) {
    Rtatami::BoundNumericPointer input(raw_input);
    const auto& shared = input->ptr;

    if (shared->sparse()) {
        auto frag = tatami::retrieve_fragmented_sparse_contents<double, int>(shared.get(), false, threads);
        const auto& store_v = frag.value;
        const auto& store_i = frag.index;

        size_t primary = shared->ncol();
        Rcpp::IntegerVector output_p(primary + 1);
        for (size_t p = 0; p < primary; ++p) {
            output_p[p + 1] = output_p[p] + store_v[p].size();
        }

        Rcpp::NumericVector output_v(output_p[primary]);
        Rcpp::IntegerVector output_i(output_p[primary]);
        size_t offset = 0;
        for (size_t p = 0; p < primary; ++p) {
            std::copy(store_v[p].begin(), store_v[p].end(), output_v.begin() + offset);
            std::copy(store_i[p].begin(), store_i[p].end(), output_i.begin() + offset);
            offset += store_v[p].size();
        }

        Rcpp::S4 output("dgCMatrix");
        output.slot("x") = output_v;
        output.slot("i") = output_i;
        output.slot("p") = output_p;
        output.slot("Dim") = Rcpp::IntegerVector::create(shared->nrow(), shared->ncol());
        return output;

    } else {
        Rcpp::NumericMatrix output(shared->nrow(), shared->ncol());
        tatami::convert_to_dense(shared.get(), false, static_cast<double*>(output.begin()), threads);
        return output;
    }
}
