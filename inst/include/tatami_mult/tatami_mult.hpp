#ifndef TATAMI_MULT_HPP
#define TATAMI_MULT_HPP

#include "tatami/tatami.hpp"

#include "dense_row.hpp"
#include "sparse_row.hpp"
#include "dense_column.hpp"
#include "sparse_column.hpp"

/**
 * @file tatami_mult.hpp
 * @brief Multiplication of **tatami** matrices.
 */

/**
 * @namespace tatami_mult
 * @brief Multiplication of **tatami** matrices.
 */
namespace tatami_mult {

/**
 * @brief Options for `multiply()`.
 */
struct Options {
    /**
     * Number of threads to use.
     */
    int num_threads = 1;

    /**
     * Whether to iterate over the preferred dimension of the larger matrix,
     * for the `multiply()` overload that accepts two `tatami::Matrix` objects.
     * This organizes the multiplication so that it only passes over the larger matrix once,
     * while passing through the other matrix multiple times. 
     * If false, the multiplication is performed with the supplied left and right matrices.
     */
    bool prefer_larger = true;

    /**
     * Whether to save the product as a column-major matrix in `output`,
     * for the `multiply()` overload that accepts two `tatami::Matrix` objects.
     * If false, the product is instead saved as a row-major matrix.
     */
    bool column_major_output = true;
};

/**
 * Compute the product `left * right`, storing the result in `output`.
 *
 * @tparam Value_ Numeric type of the matrix value.
 * @tparam Index_ Integer type of the matrix index.
 * @tparam Right_ Numeric type of the vector on the right hand side.
 * @tparam Output_ Numeric type of the output array.
 *
 * @param left A **tatami** matrix to be multiplied.
 * @param[in] right Pointer to an array of length `left.ncol()`,
 * containing the vector with which to multiply `left`.
 * @param[out] output Pointer to an array of length `left.nrow()`.
 * On output, this stores the product `left * right`.
 * @param opt Further options.
 */
template<typename Value_, typename Index_, typename Right_, typename Output_>
void multiply(const tatami::Matrix<Value_, Index_>& left, const Right_* right, Output_* output, const Options& opt) {
    if (left.sparse()) {
        if (left.prefer_rows()) {
            internal::sparse_row_vector(left, right, output, opt.num_threads);
        } else {
            internal::sparse_column_vector(left, right, output, opt.num_threads);
        }
    } else {
        if (left.prefer_rows()) {
            internal::dense_row_vector(left, right, output, opt.num_threads);
        } else {
            internal::dense_column_vector(left, right, output, opt.num_threads);
        }
    }
}

/**
 * Compute the product `t(left) * right`, storing the result in `output`.
 *
 * @tparam Left_ Numeric type of the vector on the left hand side.
 * @tparam Value_ Numeric type of the matrix value.
 * @tparam Index_ Integer type of the matrix index.
 * @tparam Output_ Numeric type of the output array.
 *
 * @param[in] left Pointer to an array of length `right.nrow()`,
 * containing the vector with which to multiply `right`.
 * @param right A **tatami** matrix to be multiplied. 
 * @param[out] output Pointer to an array of length `right.ncol()`.
 * On output, this stores the product `t(left) * right`.
 * @param opt Further options.
 */
template<typename Left_, typename Value_, typename Index_, typename Output_>
void multiply(const Left_* left, const tatami::Matrix<Value_, Index_>& right, Output_* output, const Options& opt) {
    auto tright = tatami::make_DelayedTranspose(tatami::wrap_shared_ptr(&right));
    if (tright->sparse()) {
        if (tright->prefer_rows()) {
            internal::sparse_row_vector(*tright, left, output, opt.num_threads);
        } else {
            internal::sparse_column_vector(*tright, left, output, opt.num_threads);
        }
    } else {
        if (tright->prefer_rows()) {
            internal::dense_row_vector(*tright, left, output, opt.num_threads);
        } else {
            internal::dense_column_vector(*tright, left, output, opt.num_threads);
        }
    }
}

/**
 * Compute the product `left * right[i]` for each array `i`, storing the result in `output[i]`.
 *
 * @tparam Value_ Numeric type of the matrix value.
 * @tparam Index_ Integer type of the matrix index.
 * @tparam Right_ Numeric type of the vectors on the right hand side.
 * @tparam Output_ Numeric type of the output array.
 *
 * @param left A **tatami** matrix to be multiplied.
 * @param[in] right Vector of pointers, each of which points to an array of length `left.ncol()`.
 * Each entry contains a vector with which to multiply `left`.
 * @param[out] output Vector of pointers, each of which points to an array of length `left.nrow()`.
 * On output, the `i`-th entry stores the product `left * right[i]`.
 * @param opt Further options.
 */
template<typename Value_, typename Index_, typename Right_, typename Output_>
void multiply(const tatami::Matrix<Value_, Index_>& left, const std::vector<Right_*>& right, const std::vector<Output_*>& output, const Options& opt) {
    if (left.sparse()) {
        if (left.prefer_rows()) {
            internal::sparse_row_vectors(left, right, output, opt.num_threads);
        } else {
            internal::sparse_column_vectors(left, right, output, opt.num_threads);
        }
    } else {
        if (left.prefer_rows()) {
            internal::dense_row_vectors(left, right, output, opt.num_threads);
        } else {
            internal::dense_column_vectors(left, right, output, opt.num_threads);
        }
    }
}

/**
 * Compute the product `t(left[i]) * right` for each array `i`, storing the result in `output[i]`.
 *
 * @tparam Left_ Numeric type of the vector on the left hand side.
 * @tparam Value_ Numeric type of the matrix value.
 * @tparam Index_ Integer type of the matrix index.
 * @tparam Output_ Numeric type of the output array.
 *
 * @param[in] left Vector of pointers, each of which points to an array of length `right.nrow()`.
 * Each entry contains a vector with which to multiply `right`.
 * @param right A **tatami** matrix to be multiplied. 
 * @param[out] output Vector of pointers, each of which points to an array of length `right.ncol()`.
 * On output, the `i`-th entry stores the product `t(left[i]) * right`.
 * @param opt Further options.
 */
template<typename Left_, typename Value_, typename Index_, typename Output_>
void multiply(const std::vector<Left_*>& left, const tatami::Matrix<Value_, Index_>& right, const std::vector<Output_*>& output, const Options& opt) {
    auto tright = tatami::make_DelayedTranspose(tatami::wrap_shared_ptr(&right));
    if (tright->sparse()) {
        if (tright->prefer_rows()) {
            internal::sparse_row_vectors(*tright, left, output, opt.num_threads);
        } else {
            internal::sparse_column_vectors(*tright, left, output, opt.num_threads);
        }
    } else {
        if (tright->prefer_rows()) {
            internal::dense_row_vectors(*tright, left, output, opt.num_threads);
        } else {
            internal::dense_column_vectors(*tright, left, output, opt.num_threads);
        }
    }
}

/**
 * @cond
 */
namespace internal {

template<typename LeftValue_, typename LeftIndex_, typename RightValue_, typename RightIndex_, typename Output_>
void multiply(const tatami::Matrix<LeftValue_, LeftIndex_>& left, const tatami::Matrix<RightValue_, RightIndex_>& right, Output_* output, bool column_major_out, int num_threads) {
    size_t row_shift, col_shift;
    if (column_major_out) {
        row_shift = 1;
        col_shift = left.nrow();
    } else {
        row_shift = right.ncol();
        col_shift = 1;
    }

    if (left.sparse()) {
        if (left.prefer_rows()) {
            if (right.sparse()) {
                internal::sparse_row_tatami_sparse(left, right, output, row_shift, col_shift, num_threads);
            } else {
                internal::sparse_row_tatami_dense(left, right, output, row_shift, col_shift, num_threads);
            }
        } else {
            if (right.sparse()) {
                internal::sparse_column_tatami_sparse(left, right, output, row_shift, col_shift, num_threads);
            } else {
                internal::sparse_column_tatami_dense(left, right, output, row_shift, col_shift, num_threads);
            }
        }
    } else {
        if (left.prefer_rows()) {
            if (right.sparse()) {
                internal::dense_row_tatami_sparse(left, right, output, row_shift, col_shift, num_threads);
            } else {
                internal::dense_row_tatami_dense(left, right, output, row_shift, col_shift, num_threads);
            }
        } else {
            if (right.sparse()) {
                internal::dense_column_tatami_sparse(left, right, output, row_shift, col_shift, num_threads);
            } else {
                internal::dense_column_tatami_dense(left, right, output, row_shift, col_shift, num_threads);
            }
        }
    }
}

}
/**
 * @endcond
 */

/**
 * Compute the product `left * right`, storing the result in `output`.
 *
 * @tparam RightValue_ Numeric type of the right matrix value.
 * @tparam RightIndex_ Integer type of the right matrix index.
 * @tparam LeftValue_ Numeric type of the left matrix value.
 * @tparam LeftIndex_ Integer type of the left matrix index.
 * @tparam Output_ Numeric type of the output array.
 *
 * @param left A **tatami** matrix to use in the multiplication.
 * @param right A **tatami** matrix to use in the multiplication.
 * `right.nrow()` and `left.ncol()` should be equal.
 * @param[out] output Pointer to an array of length equal to `left.nrow() * right.ncol()`. 
 * On output, this stores the matrix product as a column- or row-major matrix (see `Options::column_major_output`).
 * @param opt Further options.
 */
template<typename LeftValue_, typename LeftIndex_, typename RightValue_, typename RightIndex_, typename Output_>
void multiply(const tatami::Matrix<LeftValue_, LeftIndex_>& left, const tatami::Matrix<RightValue_, RightIndex_>& right, Output_* output, const Options& opt) {
    if (opt.prefer_larger) {
        if (left.nrow() < right.ncol()) {
            auto tright = tatami::make_DelayedTranspose(tatami::wrap_shared_ptr(&right));
            auto tleft = tatami::make_DelayedTranspose(tatami::wrap_shared_ptr(&left));
            internal::multiply(*tright, *tleft, output, !opt.column_major_output, opt.num_threads);
            return;
        }
    }

    internal::multiply(left, right, output, opt.column_major_output, opt.num_threads);
}

}

#endif
