#include "Rtatami.h"
#include "Rcpp.h"
#include "tatami_stats/tatami_stats.hpp"
#include "tatami_mult/tatami_mult.hpp"

#include <vector>
#include <cstddef>
#include <stdexcept>

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

    auto wrk = tatami::new_extractor<false, false>(*(input->ptr), row, false);
    Rcpp::NumericVector output(row ? NC : NR);
    auto optr = static_cast<double*>(output.begin());
    auto ptr = wrk->fetch(i - 1, optr);
    tatami::copy_n(ptr, output.size(), optr);
    return output;
}

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector tatami_sums(SEXP raw_input, bool row, int threads) {
    tatami_stats::sums::Options opt;
    opt.num_threads = threads;
    if (threads < 1) {
        throw std::runtime_error("'threads' should be a positive integer");
    }

    Rtatami::BoundNumericPointer input(raw_input);
    const auto NR = input->ptr->nrow();
    const auto NC = input->ptr->ncol();

    Rcpp::NumericVector output(row ? NR : NC);
    tatami_stats::sums::apply(row, *(input->ptr), static_cast<double*>(output.begin()), opt);
    return output;
}

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector tatami_sums_by_group(SEXP raw_input, Rcpp::IntegerVector group, bool row, int threads) {
    tatami_stats::grouped_sums::Options opt;
    opt.num_threads = threads;
    if (threads < 1) {
        throw std::runtime_error("'threads' should be a positive integer");
    }

    Rtatami::BoundNumericPointer input(raw_input);
    const auto NR = input->ptr->nrow();
    const auto NC = input->ptr->ncol();

    std::vector<int> group_m1(group.begin(), group.end());
    for (auto& x : group_m1) {
        if (x <= 0) {
            throw std::runtime_error("'group' should contain positive group identifiers");
        }
        --x;
    }
    const auto total_groups = tatami_stats::total_groups(group_m1.data(), group_m1.size());

    const auto expected = (row ? NC : NR);
    if (!sanisizer::is_equal(expected, group_m1.size())) {
        throw std::runtime_error("'group' should have length equal to the appropriate dimension extent");
    }

    const auto stride = (row ? NR : NC);
    Rcpp::NumericMatrix output(stride, total_groups);
    auto ptrs = sanisizer::create<std::vector<double*> >(total_groups);
    for (auto i = static_cast<decltype(total_groups)>(0); i < total_groups; ++i) {
        ptrs[i] = output.begin() + sanisizer::product_unsafe<std::size_t>(i, stride);
    }
    tatami_stats::grouped_sums::apply(row, *(input->ptr), group_m1.data(), total_groups, ptrs.data(), opt);

    if (row) {
        return output;
    } else {
        return Rcpp::transpose(output);
    }
}

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector tatami_medians(SEXP raw_input, bool row, int threads) {
    tatami_stats::medians::Options opt;
    opt.num_threads = threads;
    if (threads < 1) {
        throw std::runtime_error("'threads' should be a positive integer");
    }

    Rtatami::BoundNumericPointer input(raw_input);
    const auto NR = input->ptr->nrow();
    const auto NC = input->ptr->ncol();

    Rcpp::NumericVector output(row ? NR : NC);
    tatami_stats::medians::apply(row, *(input->ptr), static_cast<double*>(output.begin()), opt);
    return output;
}

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector tatami_nan_counts(SEXP raw_input, bool row, int threads) {
    tatami_stats::counts::nan::Options opt;
    opt.num_threads = threads;
    if (threads < 1) {
        throw std::runtime_error("'threads' should be a positive integer");
    }

    Rtatami::BoundNumericPointer input(raw_input);
    const auto NR = input->ptr->nrow();
    const auto NC = input->ptr->ncol();

    Rcpp::NumericVector output(row ? NR : NC);
    tatami_stats::counts::nan::apply(row, *(input->ptr), static_cast<double*>(output.begin()), opt);
    return output;
}

//[[Rcpp::export(rng=false)]]
SEXP tatami_realize(SEXP raw_input, int threads) {
    if (threads < 1) {
        throw std::runtime_error("'threads' should be a positive integer");
    }

    Rtatami::BoundNumericPointer input(raw_input);
    const auto& shared = input->ptr;

    if (shared->sparse()) {
        auto frag = tatami::retrieve_fragmented_sparse_contents<double, int>(shared.get(), false, threads);
        const auto& store_v = frag.value;
        const auto& store_i = frag.index;

        auto primary = shared->ncol();
        Rcpp::IntegerVector output_p(sanisizer::sum<decltype(std::declval<Rcpp::IntegerVector>().size())>(primary, 1));
        for (decltype(primary) p = 0; p < primary; ++p) {
            output_p[p + 1] = output_p[p] + store_v[p].size();
        }

        auto last_p = output_p[primary];
        auto output_v = sanisizer::create<Rcpp::NumericVector>(last_p);
        auto output_i = sanisizer::create<Rcpp::IntegerVector>(last_p);
        decltype(last_p) offset = 0;
        for (decltype(primary) p = 0; p < primary; ++p) {
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
Rcpp::NumericVector tatami_multiply_vector(SEXP raw_input, Rcpp::NumericVector other, bool right, int threads) {
    if (threads < 1) {
        throw std::runtime_error("'threads' should be a positive integer");
    }

    Rtatami::BoundNumericPointer input(raw_input);
    const auto& shared = input->ptr;

    tatami_mult::Options opt;
    opt.num_threads = threads; 
    if (right) {
        if (!sanisizer::is_equal(other.size(), shared->ncol())) {
            throw std::runtime_error("length of vector does not match the number of columns of 'x'");
        }
        Rcpp::NumericVector output(shared->nrow());
        tatami_mult::multiply(*shared, static_cast<const double*>(other.begin()), static_cast<double*>(output.begin()), opt);
        return output;
    } else {
        if (!sanisizer::is_equal(other.size(), shared->nrow())) {
            throw std::runtime_error("length of vector does not match the number of rows of 'x'");
        }
        Rcpp::NumericVector output(shared->ncol());
        tatami_mult::multiply(static_cast<const double*>(other.begin()), *shared, static_cast<double*>(output.begin()), opt);
        return output;
    }
}

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector tatami_multiply_columns(SEXP raw_input, Rcpp::NumericMatrix other, bool right, int threads) {
    if (threads < 1) {
        throw std::runtime_error("'threads' should be a positive integer");
    }

    Rtatami::BoundNumericPointer input(raw_input);
    const auto& shared = input->ptr;

    tatami_mult::Options opt;
    opt.num_threads = threads; 

    if (right) {
        auto common_dim = other.rows();
        if (!sanisizer::is_equal(common_dim, shared->ncol())) {
            throw std::runtime_error("rows of 'vals' does not match the number of columns of 'x'");
        }

        auto num_other = other.cols();
        auto out_nrow = shared->nrow();
        sanisizer::product<std::size_t>(out_nrow, num_other); // check outside the loop so that we can do unsafe products inside the loop.
        Rcpp::NumericMatrix output(sanisizer::cast<decltype(num_other)>(out_nrow), num_other);

        auto other_ptrs = sanisizer::create<std::vector<const double*> >(num_other);
        auto out_ptrs = sanisizer::create<std::vector<double*> >(num_other);
        double* outptr = output.begin();
        const double* otherptr = other.begin();
        for (decltype(num_other) o = 0; o < num_other; ++o) {
            out_ptrs[o] = outptr + sanisizer::product_unsafe<std::size_t>(o, out_nrow);
            other_ptrs[o] = otherptr + sanisizer::product_unsafe<std::size_t>(o, common_dim);
        }

        tatami_mult::multiply(*shared, other_ptrs, out_ptrs, opt);
        return output;

    } else {
        // In this case, we assume that 'other' is in its transposed form so that we get
        // column-major access to what was previously the rows in the original matrix.
        auto common_dim = other.rows();
        if (!sanisizer::is_equal(common_dim, shared->nrow())) {
            throw std::runtime_error("columns of 'vals' does not match the number of rows of 'x'");
        }

        auto num_other = other.cols();
        auto out_nrow = shared->ncol();
        sanisizer::product<std::size_t>(out_nrow, num_other); // check outside the loop so that we can do unsafe products inside the loop.
        Rcpp::NumericMatrix output(sanisizer::cast<decltype(num_other)>(out_nrow), num_other); // also transposed, so that we can write conveniently to a column-major format.

        auto other_ptrs = sanisizer::create<std::vector<const double*> >(num_other);
        auto out_ptrs = sanisizer::create<std::vector<double*> >(num_other);
        double* outptr = output.begin();
        const double* otherptr = other.begin();
        for (decltype(num_other) o = 0; o < num_other; ++o) {
            out_ptrs[o] = outptr + sanisizer::product_unsafe<std::size_t>(o, out_nrow);
            other_ptrs[o] = otherptr + sanisizer::product_unsafe<std::size_t>(o, common_dim);
        }

        tatami_mult::multiply(other_ptrs, *shared, out_ptrs, opt);
        return output;
    } 
}

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector tatami_multiply_matrix(SEXP raw_input, SEXP more_input, bool right, int threads) {
    if (threads < 1) {
        throw std::runtime_error("'threads' should be a positive integer");
    }

    Rtatami::BoundNumericPointer input(raw_input);
    Rtatami::BoundNumericPointer input2(more_input);

    const auto& shared = (right ? input->ptr : input2->ptr);
    const auto& shared2 = (right ? input2->ptr : input->ptr);

    if (shared->ncol() != shared2->nrow()) {
        throw std::runtime_error("inconsistent common dimensions for matrix multiplication");
    }

    tatami_mult::Options opt;
    opt.num_threads = threads; 
    Rcpp::NumericMatrix output(shared->nrow(), shared2->ncol());
    tatami_mult::multiply(*shared, *shared2, static_cast<double*>(output.begin()), opt);

    return output;
}
