#include "Rtatami.h"
#include "Rcpp.h"

#include "tatami/tatami.hpp"
#include "tatami_mult/tatami_mult.hpp"
#include "sanisizer/sanisizer.hpp"

#include <vector>
#include <cstddef>
#include <stdexcept>

#include "utils.h"

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
        auto output = tatami::create_container_of_Index_size<Rcpp::NumericVector>(shared->nrow());
        tatami_mult::multiply_with_single_vector(*shared, static_cast<const double*>(other.begin()), static_cast<double*>(output.begin()), opt);
        return output;

    } else {
        if (!sanisizer::is_equal(other.size(), shared->nrow())) {
            throw std::runtime_error("length of vector does not match the number of rows of 'x'");
        }
        auto output = tatami::create_container_of_Index_size<Rcpp::NumericVector>(shared->ncol());
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
        auto output = create_output_matrix(out_nrow, num_other);
        const auto optr = static_cast<double*>(output.begin()); 

        const auto rptr = static_cast<const double*>(other.begin());
        sanisizer::product<std::size_t>(num_other, common_dim); // check outside the functor so that we can do unsafe products inside the functor.
        auto fetch_other_col = [&](std::size_t rcol) -> const double* {
            return rptr + sanisizer::product_unsafe<std::size_t>(rcol, common_dim);
        };

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
        auto output = create_output_matrix(num_other, out_nrow);
        const auto optr = static_cast<double*>(output.begin()); 

        const auto rptr = static_cast<const double*>(other.begin());
        sanisizer::product<std::size_t>(num_other, common_dim); // check outside the functor so that we can do unsafe products inside the functor.
        auto fetch_other_row = [&](std::size_t rrow) -> const double* {
            return rptr + sanisizer::product_unsafe<std::size_t>(rrow, num_other);
        };

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

