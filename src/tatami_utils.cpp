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
    tatami_stats::SumOptions opt;
    opt.num_threads = threads;
    if (threads < 1) {
        throw std::runtime_error("'threads' should be a positive integer");
    }

    Rtatami::BoundNumericPointer input(raw_input);
    const auto NR = input->ptr->nrow();
    const auto NC = input->ptr->ncol();

    Rcpp::NumericVector output(row ? NR : NC);
    tatami_stats::sum(row, *(input->ptr), static_cast<double*>(output.begin()), opt);
    return output;
}

//[[Rcpp::export(rng=false)]]
Rcpp::List tatami_variances(SEXP raw_input, bool row, int threads) {
    tatami_stats::VarianceOptions opt;
    opt.num_threads = threads;
    if (threads < 1) {
        throw std::runtime_error("'threads' should be a positive integer");
    }

    Rtatami::BoundNumericPointer input(raw_input);
    const auto NR = input->ptr->nrow();
    const auto NC = input->ptr->ncol();

    Rcpp::NumericVector output_means(row ? NR : NC);
    Rcpp::NumericVector output_variances(row ? NR : NC);

    tatami_stats::VarianceBuffers<double> buffers;
    buffers.mean = output_means.begin();
    buffers.variance = output_variances.begin();
    tatami_stats::variance(row, *(input->ptr), buffers, opt);

    return Rcpp::List::create(
        Rcpp::Named("mean") = output_means,
        Rcpp::Named("variance") = output_variances
    );
}

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector tatami_medians(SEXP raw_input, bool row, int threads) {
    tatami_stats::MedianOptions opt;
    opt.num_threads = threads;
    if (threads < 1) {
        throw std::runtime_error("'threads' should be a positive integer");
    }

    Rtatami::BoundNumericPointer input(raw_input);
    const auto NR = input->ptr->nrow();
    const auto NC = input->ptr->ncol();

    Rcpp::NumericVector output(row ? NR : NC);
    tatami_stats::median(row, *(input->ptr), static_cast<double*>(output.begin()), opt);
    return output;
}

//[[Rcpp::export(rng=false)]]
Rcpp::List tatami_ranges(SEXP raw_input, bool row, int threads) {
    tatami_stats::RangeOptions opt;
    opt.num_threads = threads;
    if (threads < 1) {
        throw std::runtime_error("'threads' should be a positive integer");
    }

    Rtatami::BoundNumericPointer input(raw_input);
    const auto NR = input->ptr->nrow();
    const auto NC = input->ptr->ncol();

    const auto expected = row ? NR : NC;
    Rcpp::NumericVector minout(expected);
    Rcpp::NumericVector maxout(expected);

    tatami_stats::RangeBuffers<double> buf;
    buf.minimum = minout.begin();
    buf.maximum = maxout.begin();

    tatami_stats::range(row, *(input->ptr), buf, opt);
    return Rcpp::List::create(
        Rcpp::Named("minimum") = minout,
        Rcpp::Named("maximum") = maxout
    );
}

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector tatami_nan_counts(SEXP raw_input, bool row, int threads) {
    tatami_stats::CountOptions opt;
    opt.num_threads = threads;
    if (threads < 1) {
        throw std::runtime_error("'threads' should be a positive integer");
    }

    Rtatami::BoundNumericPointer input(raw_input);
    const auto NR = input->ptr->nrow();
    const auto NC = input->ptr->ncol();

    Rcpp::NumericVector output(row ? NR : NC);
    tatami_stats::count(
        row,
        *(input->ptr),
        static_cast<double*>(output.begin()),
        [](const double val) -> bool { return std::isnan(val); },
        opt
    );
    return output;
}

static std::vector<int> extract_groups(const Rcpp::IntegerVector& group, int num_groups) {
    std::vector<int> group_m1(group.begin(), group.end());
    for (auto& x : group_m1) {
        if (x <= 0) {
            throw std::runtime_error("'group' should contain positive group identifiers");
        } else if (x > num_groups) {
            throw std::runtime_error("entries of 'group' should be less than 'num.groups'");
        }
        --x;
    }
    return group_m1;
}

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector tatami_sums_by_group(SEXP raw_input, Rcpp::IntegerVector group, int num_groups, bool row, int threads) {
    tatami_stats::GroupSumOptions opt;
    opt.num_threads = threads;
    if (threads < 1) {
        throw std::runtime_error("'threads' should be a positive integer");
    }

    Rtatami::BoundNumericPointer input(raw_input);
    const auto NR = input->ptr->nrow();
    const auto NC = input->ptr->ncol();

    auto group_m1 = extract_groups(group, num_groups);
    const auto expected = (row ? NC : NR);
    if (!sanisizer::is_equal(expected, group_m1.size())) {
        throw std::runtime_error("'group' should have length equal to the appropriate dimension extent");
    }

    const auto stride = (row ? NR : NC);
    Rcpp::NumericMatrix output(stride, num_groups);
    auto ptrs = sanisizer::create<std::vector<double*> >(num_groups);
    for (auto i = static_cast<decltype(num_groups)>(0); i < num_groups; ++i) {
        ptrs[i] = output.begin() + sanisizer::product_unsafe<std::size_t>(i, stride);
    }

    tatami_stats::group_sum(row, *(input->ptr), group_m1.data(), num_groups, ptrs, opt);

    if (row) {
        return output;
    } else {
        return Rcpp::transpose(output);
    }
}

//[[Rcpp::export(rng=false)]]
Rcpp::List tatami_variances_by_group(SEXP raw_input, Rcpp::IntegerVector group, int num_groups, bool row, int threads) {
    tatami_stats::GroupVarianceOptions opt;
    opt.num_threads = threads;
    if (threads < 1) {
        throw std::runtime_error("'threads' should be a positive integer");
    }

    Rtatami::BoundNumericPointer input(raw_input);
    const auto NR = input->ptr->nrow();
    const auto NC = input->ptr->ncol();

    auto group_m1 = extract_groups(group, num_groups);
    const auto expected = (row ? NC : NR);
    if (!sanisizer::is_equal(expected, group_m1.size())) {
        throw std::runtime_error("'group' should have length equal to the appropriate dimension extent");
    }

    const auto stride = (row ? NR : NC);
    Rcpp::NumericMatrix moutput(stride, num_groups), voutput(stride, num_groups);
    tatami_stats::GroupVarianceBuffers<double> buf;
    buf.mean.reserve(num_groups);
    buf.variance.reserve(num_groups);
    sanisizer::product<std::size_t>(stride, num_groups); // check that all offset calculations are safe.
    for (auto i = static_cast<decltype(num_groups)>(0); i < num_groups; ++i) {
        const auto offset = sanisizer::product_unsafe<std::size_t>(i, stride);
        buf.mean.push_back(moutput.begin() + offset);
        buf.variance.push_back(voutput.begin() + offset);
    }

    tatami_stats::group_variance(row, *(input->ptr), group_m1.data(), num_groups, buf, opt);

    if (row) {
        return Rcpp::List::create(
            Rcpp::Named("mean") = moutput,
            Rcpp::Named("variance") = voutput 
        );
    } else {
        return Rcpp::List::create(
            Rcpp::Named("mean") = Rcpp::transpose(moutput),
            Rcpp::Named("variance") = Rcpp::transpose(voutput)
        );
    }
}

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector tatami_medians_by_group(SEXP raw_input, Rcpp::IntegerVector group, int num_groups, bool row, int threads) {
    tatami_stats::GroupMedianOptions opt;
    opt.num_threads = threads;
    if (threads < 1) {
        throw std::runtime_error("'threads' should be a positive integer");
    }

    Rtatami::BoundNumericPointer input(raw_input);
    const auto NR = input->ptr->nrow();
    const auto NC = input->ptr->ncol();

    auto group_m1 = extract_groups(group, num_groups);
    const auto expected = (row ? NC : NR);
    if (!sanisizer::is_equal(expected, group_m1.size())) {
        throw std::runtime_error("'group' should have length equal to the appropriate dimension extent");
    }

    const auto stride = (row ? NR : NC);
    Rcpp::NumericMatrix output(stride, num_groups);
    auto ptrs = sanisizer::create<std::vector<double*> >(num_groups);
    sanisizer::product<std::size_t>(stride, num_groups); // check that product inside the loop is safe.
    for (auto i = static_cast<decltype(num_groups)>(0); i < num_groups; ++i) {
        ptrs[i] = output.begin() + sanisizer::product_unsafe<std::size_t>(i, stride);
    }

    tatami_stats::group_median(row, *(input->ptr), group_m1.data(), num_groups, ptrs, opt);

    if (row) {
        return output;
    } else {
        return Rcpp::transpose(output);
    }
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
    tatami_mult::MultiplyWithSingleVectorOptions opt;
    tatami_mult::set_num_threads(opt, threads); 

    Rtatami::BoundNumericPointer input(raw_input);
    const auto& shared = input->ptr;

    if (right) {
        if (!sanisizer::is_equal(other.size(), shared->ncol())) {
            throw std::runtime_error("length of vector does not match the number of columns of 'x'");
        }
        Rcpp::NumericVector output(shared->nrow());
        tatami_mult::multiply_with_single_vector(*shared, static_cast<const double*>(other.begin()), static_cast<double*>(output.begin()), opt);
        return output;

    } else {
        if (!sanisizer::is_equal(other.size(), shared->nrow())) {
            throw std::runtime_error("length of vector does not match the number of rows of 'x'");
        }
        Rcpp::NumericVector output(shared->ncol());
        tatami_mult::multiply_with_single_vector(static_cast<const double*>(other.begin()), *shared, static_cast<double*>(output.begin()), opt);
        return output;
    }
}

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector tatami_multiply_columns(SEXP raw_input, Rcpp::NumericMatrix other, bool right, int threads) {
    if (threads < 1) {
        throw std::runtime_error("'threads' should be a positive integer");
    }
    Rtatami::BoundNumericPointer input(raw_input);

    if (right) {
        const auto common_dim = other.rows();
        const auto num_other = other.cols();

        const auto& mat = *(input->ptr);
        if (!sanisizer::is_equal(common_dim, mat.ncol())) {
            throw std::runtime_error("number of rows of 'vals' does not match the number of columns of 'x'");
        }

        const auto out_nrow = mat.nrow();
        Rcpp::NumericMatrix output(out_nrow, num_other);
        const auto optr = static_cast<double*>(output.begin()); 

        const auto rptr = static_cast<const double*>(other.begin());
        sanisizer::product<std::size_t>(num_other, common_dim); // check outside the functor so that we can do unsafe products inside the functor.
        auto fetch_other_col = [&](std::size_t rcol) -> const double* { return rptr + sanisizer::product_unsafe<std::size_t>(rcol, common_dim); };

        const bool is_sparse = mat.is_sparse(); 
        if (mat.prefer_rows()) {
            if (is_sparse){
                tatami_mult::MultiplySparseRowWithDenseColumnMatrixToColumnOutputOptions opt;
                opt.num_threads = threads;
                tatami_mult::multiply_sparse_row_with_dense_column_matrix_to_column_output(mat, num_other, fetch_other_col, optr, opt);
            } else {
                tatami_mult::MultiplyDenseRowWithDenseColumnMatrixToColumnOutputOptions opt;
                opt.num_threads = threads;
                tatami_mult::multiply_dense_row_with_dense_column_matrix_to_column_output(mat, num_other, fetch_other_col, optr, opt);
            }

        } else {
            if (is_sparse){
                tatami_mult::MultiplySparseColumnWithDenseColumnMatrixToColumnOutputOptions opt;
                opt.num_threads = threads;
                tatami_mult::multiply_sparse_column_with_dense_column_matrix_to_column_output(mat, num_other, fetch_other_col, optr, opt);
            } else {
                tatami_mult::MultiplyDenseColumnWithDenseColumnMatrixToColumnOutputOptions opt;
                opt.num_threads = threads;
                tatami_mult::multiply_dense_column_with_dense_column_matrix_to_column_output(mat, num_other, fetch_other_col, optr, opt);
            }
        }

        return output;

    } else {
        // Effectively transposing 'other' here so that it can lie on the right as a row-major matrix.
        const auto common_dim = other.cols();
        const auto num_other = other.rows();

        tatami::DelayedTranspose<double, int> alt(input->ptr);
        if (!sanisizer::is_equal(common_dim, alt.ncol())) {
            throw std::runtime_error("number of columns of 'vals' does not match the number of rows of 'x'");
        }

        // We'll treat the output as row-major in our (conceptually) everything-is-transposed state.
        const auto out_nrow = alt.nrow();
        Rcpp::NumericMatrix output(num_other, out_nrow);
        const auto optr = static_cast<double*>(output.begin()); 

        const auto rptr = static_cast<const double*>(other.begin());
        sanisizer::product<std::size_t>(num_other, common_dim); // check outside the functor so that we can do unsafe products inside the functor.
        auto fetch_other_row = [&](std::size_t rrow) -> const double* { return rptr + sanisizer::product_unsafe<std::size_t>(rrow, num_other); };

        const bool is_sparse = alt.is_sparse(); 
        if (alt.prefer_rows()) {
            if (is_sparse){
                tatami_mult::MultiplySparseRowWithDenseRowMatrixToRowOutputOptions opt;
                opt.num_threads = threads;
                tatami_mult::multiply_sparse_row_with_dense_row_matrix_to_row_output(alt, num_other, fetch_other_row, optr, opt);
            } else {
                tatami_mult::MultiplyDenseRowWithDenseRowMatrixToRowOutputOptions opt;
                opt.num_threads = threads;
                tatami_mult::multiply_dense_row_with_dense_row_matrix_to_row_output(alt, num_other, fetch_other_row, optr, opt);
            }

        } else {
            if (is_sparse){
                tatami_mult::MultiplySparseColumnWithDenseRowMatrixToRowOutputOptions opt;
                opt.num_threads = threads;
                tatami_mult::multiply_sparse_column_with_dense_row_matrix_to_row_output(alt, num_other, fetch_other_row, optr, opt);
            } else {
                tatami_mult::MultiplyDenseColumnWithDenseRowMatrixToRowOutputOptions opt;
                opt.num_threads = threads;
                tatami_mult::multiply_dense_column_with_dense_row_matrix_to_row_output(alt, num_other, fetch_other_row, optr, opt);
            }
        }

        return output;
    }
}

//[[Rcpp::export(rng=false)]]
Rcpp::NumericMatrix tatami_multiply_matrix(SEXP raw_input1, SEXP raw_input2, bool right, int threads) {
    if (threads < 1) {
        throw std::runtime_error("'threads' should be a positive integer");
    }
    tatami_mult::MultiplyWithMatrixOptions opt;
    tatami_mult::set_num_threads(opt, threads); 

    Rtatami::BoundNumericPointer input1(raw_input1);
    Rtatami::BoundNumericPointer input2(raw_input2);
    const auto& lmat = (right ? input1->ptr : input2->ptr);
    const auto& rmat = (right ? input2->ptr : input1->ptr);
    if (lmat->ncol() != rmat->nrow()) {
        throw std::runtime_error("inconsistent common dimensions for matrix multiplication");
    }

    Rcpp::NumericMatrix output(lmat->nrow(), rmat->ncol());
    sanisizer::product<std::size_t>(lmat->nrow(), rmat->ncol()); // check that the contents of the array are accessible via a size_t.
    tatami_mult::multiply_with_matrix(*lmat, *rmat, static_cast<double*>(output.begin()), false, opt);
    return output;
}
