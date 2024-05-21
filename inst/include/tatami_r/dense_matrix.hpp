#ifndef TATAMI_R_DENSE_MATRIX_HPP
#define TATAMI_R_DENSE_MATRIX_HPP

#include "tatami/tatami.hpp"
#include <algorithm>

namespace tatami_r { 

template<typename InputValue_, class InputObject_,  typename CachedValue_>
void parse_dense_matrix_internal(const InputObject_& y, bool row, CachedValue_* cache, size_t start_row, size_t start_col, size_t num_rows, size_t num_cols) {
    auto input = static_cast<const InputValue_*>(y.begin()) + start_row + start_col * static_cast<size_t>(y.rows());

    if (row) {
        // y is a column-major matrix, but transpose() expects a row-major
        // input, so we just conceptually transpose it.
        tatami::transpose(input, num_cols, num_rows, y.rows(), cache, num_cols);
    } else {
        for (size_t c = 0; c < num_cols; ++c) {
            std::copy_n(input, num_rows, cache);
            input += y.rows();
            cache += num_rows;
        }
    }
}

template<typename CachedValue_>
void parse_dense_matrix(const Rcpp::RObject& seed, bool row, CachedValue_* cache, size_t start_row, size_t start_col, size_t num_rows, size_t num_cols) {
    auto stype = seed.sexp_type();
    if (stype == REALSXP) {
        Rcpp::NumericMatrix y(seed);
        parse_dense_matrix_internal<double>(y, row, cache, start_row, start_col, num_rows, num_cols);
    } else if (stype == INTSXP) {
        Rcpp::IntegerMatrix y(seed);
        parse_dense_matrix_internal<int>(y, row, cache, start_row, start_col, num_rows, num_cols);
    } else if (stype == LGLSXP) {
        Rcpp::LogicalMatrix y(seed);
        parse_dense_matrix_internal<int>(y, row, cache, start_row, start_col, num_rows, num_cols);
    } else {
        throw std::runtime_error("unsupported SEXP type (" + std::to_string(stype) + ") from the matrix returned by 'extract_array'");
    }
}

}

#endif
