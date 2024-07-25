#include "Rtatami.h"
#include "Rcpp.h"
#include "tatami_stats/tatami_stats.hpp"
#include "tatami_mult/tatami_mult.hpp"

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

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector tatami_multiply_vector(SEXP raw_input, Rcpp::NumericVector other, bool right, int num_threads) {
    Rtatami::BoundNumericPointer input(raw_input);
    const auto& shared = input->ptr;

    tatami_mult::Options opt;
    opt.num_threads = num_threads; 
    if (right) {
        if (static_cast<int>(other.size()) != shared->ncol()) {
            throw std::runtime_error("length of vector does not match the number of columns of 'x'");
        }
        Rcpp::NumericVector output(shared->nrow());
        tatami_mult::multiply(*shared, static_cast<const double*>(other.begin()), static_cast<double*>(output.begin()), opt);
        return output;
    } else {
        if (static_cast<int>(other.size()) != shared->nrow()) {
            throw std::runtime_error("length of vector does not match the number of rows of 'x'");
        }
        Rcpp::NumericVector output(shared->ncol());
        tatami_mult::multiply(static_cast<const double*>(other.begin()), *shared, static_cast<double*>(output.begin()), opt);
        return output;
    }
}

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector tatami_multiply_columns(SEXP raw_input, Rcpp::NumericMatrix other, bool right, int num_threads) {
    Rtatami::BoundNumericPointer input(raw_input);
    const auto& shared = input->ptr;

    tatami_mult::Options opt;
    opt.num_threads = num_threads; 

    if (right) {
        size_t common_dim = other.rows();
        if (common_dim != static_cast<size_t>(shared->ncol())) {
            throw std::runtime_error("rows of 'vals' does not match the number of columns of 'x'");
        }

        size_t num_other = other.cols();
        size_t out_nrow = shared->nrow();
        Rcpp::NumericMatrix output(out_nrow, num_other);
        std::vector<double*> other_ptrs(num_other), out_ptrs(num_other);
        size_t other_offset = 0, out_offset = 0;
        for (size_t o = 0; o < num_other; ++o, other_offset += common_dim, out_offset += out_nrow) {
            out_ptrs[o] = output.begin() + out_offset;
            other_ptrs[o] = other.begin() + other_offset;
        }

        tatami_mult::multiply(*shared, other_ptrs, out_ptrs, opt);
        return output;

    } else {
        // In this case, we assume that 'other' is in its transposed form so that we get
        // column-major access to what was previously the rows in the original matrix.
        size_t common_dim = other.rows();
        if (common_dim != static_cast<size_t>(shared->nrow())) {
            throw std::runtime_error("columns of 'vals' does not match the number of rows of 'x'");
        }

        size_t num_other = other.cols();
        size_t out_nrow = shared->ncol();
        Rcpp::NumericMatrix output(out_nrow, num_other); // also transposed, so that we can write conveniently to a column-major format.
        std::vector<double*> other_ptrs(num_other), out_ptrs(num_other);
        size_t other_offset = 0, out_offset = 0;
        for (size_t o = 0; o < num_other; ++o, other_offset += common_dim, out_offset += out_nrow) {
            out_ptrs[o] = output.begin() + out_offset;
            other_ptrs[o] = other.begin() + other_offset;
        }

        tatami_mult::multiply(other_ptrs, *shared, out_ptrs, opt);
        return output;
    } 
}

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector tatami_multiply_matrix(SEXP raw_input, SEXP more_input, bool right, int num_threads) {
    Rtatami::BoundNumericPointer input(raw_input);
    Rtatami::BoundNumericPointer input2(more_input);

    const auto& shared = (right ? input->ptr : input2->ptr);
    const auto& shared2 = (right ? input2->ptr : input->ptr);

    if (shared->ncol() != shared2->nrow()) {
        throw std::runtime_error("inconsistent common dimensions for matrix multiplication");
    }

    tatami_mult::Options opt;
    opt.num_threads = num_threads; 
    Rcpp::NumericMatrix output(shared->nrow(), shared2->ncol());
    tatami_mult::multiply(*shared, *shared2, static_cast<double*>(output.begin()), opt);

    return output;
}
