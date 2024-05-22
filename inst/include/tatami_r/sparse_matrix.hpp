#ifndef TATAMI_R_SPARSE_MATRIX_HPP
#define TATAMI_R_SPARSE_MATRIX_HPP

#include "utils.hpp"
#include "tatami/tatami.hpp"
#include <type_traits>

/**
 * @file sparse_matrix.hpp
 * @brief Parse sparse matrices from block processing.
 */

namespace tatami_r { 

/**
 * Parse the contents of a `SVT_SparseMatrix` from the **DelayedArray** package.
 * This accounts for different versions of the class definition, different types of the values, and the presence of lacunar leaf nodes.
 *
 * @tparam Function_ Function to be applied at each leaf node.
 *
 * @param matrix The `SVT_SparseMatrix` object.
 * @param fun Function to apply to each leaf node, accepting four arguments:
 * .
 * 1. `c`, an integer specifying the index of the leaf node, i.e., the column index. 
 * 2. `indices`, an `Rcpp::IntegerVector` containing the sorted, zero-based indices of the structural non-zero elements in this node (i.e., column).
 * 3. `all_ones`, a boolean indicating whether all values in this node/column are equal to 1.
 * 4. `values`, an `Rcpp::IntegerVector`, `Rcpp::LogicalVector` or `Rcpp::NumericVector` containing the values of the structural non-zeros.
 *    This should be of the same length as `indices`.
 *    It should be ignored if `all_ones = true`.
 * .
 * The return value of this function is ignored.
 * Note that `fun` may not be called for all `c` - if leaf nodes do not contain any data, they will be skipped.
 */
template<class Function_>
void parse_SVT_SparseMatrix(Rcpp::RObject matrix, Function_ fun) {
    Rcpp::RObject raw_svt = matrix.slot("SVT");
    if (raw_svt == R_NilValue) {
        return;
    }

    Rcpp::IntegerVector svt_version(matrix.slot(".svt_version"));
    if (svt_version.size() != 1) {
        auto ctype = get_class_name(matrix);
        throw std::runtime_error("'.svt_version' slot of a " + ctype + " should be an integer scalar");
    }
    int version = svt_version[0];
    int index_x = (version == 0 ? 0 : 1);
    int value_x = (version == 0 ? 1 : 0);

    Rcpp::List svt(raw_svt);
    int NC = svt.size();

    for (int c = 0; c < NC; ++c) {
        Rcpp::RObject raw_inner(svt[c]);
        if (raw_inner == R_NilValue) {
            continue;
        }

        Rcpp::List inner(raw_inner);
        if (inner.size() != 2) {
            auto ctype = get_class_name(matrix);
            throw std::runtime_error("each entry of the 'SVT' slot of a " + ctype + " object should be a list of length 2 or NULL");
        }

        // Verify type to ensure that we're not making a view on a temporary array.
        Rcpp::RObject raw_indices = inner[index_x];
        if (raw_indices.sexp_type() != INTSXP) {
            auto ctype = get_class_name(matrix);
            throw std::runtime_error("indices of each element of the 'SVT' slot in a " + ctype + " object should be an integer vector");
        }
        Rcpp::IntegerVector curindices(raw_indices);
        auto nnz = curindices.size();

        Rcpp::RObject raw_values(inner[value_x]);
        auto vsexp = raw_values.sexp_type();
        bool has_values = raw_values != R_NilValue;
        Rcpp::IntegerVector curvalues_i;
        Rcpp::NumericVector curvalues_n;
        Rcpp::LogicalVector curvalues_l;

        if (has_values) {
            decltype(nnz) vsize;
            switch (vsexp) {
                case INTSXP:
                    curvalues_i = Rcpp::IntegerVector(raw_values);
                    vsize = curvalues_i.size();
                    break;
                case REALSXP:
                    curvalues_n = Rcpp::NumericVector(raw_values);
                    vsize = curvalues_n.size();
                    break;
                case LGLSXP:
                    curvalues_l = Rcpp::LogicalVector(raw_values);
                    vsize = curvalues_l.size();
                    break;
                default:
                    {
                        auto ctype = get_class_name(matrix);
                        throw std::runtime_error("value vector of an element of the 'SVT' slot in a " + ctype + " object is not a numeric or logical type");
                    }
            }

            if (nnz != vsize) {
                auto ctype = get_class_name(matrix);
                throw std::runtime_error("both vectors of an element of the 'SVT' slot in a " + ctype + " object should have the same length");
            }
        }

        switch (vsexp) {
            case INTSXP:
                fun(c, curindices, !has_values, curvalues_i);
                break;
            case REALSXP:
                fun(c, curindices, !has_values, curvalues_n);
                break;
            default:
                fun(c, curindices, !has_values, curvalues_l);
                break;
        }
    }
}

/**
 * @cond
 */
template<typename CachedValue_, typename CachedIndex_, typename Index_>
void parse_sparse_matrix(
    Rcpp::RObject matrix,
    bool row,
    std::vector<CachedValue_*>& value_ptrs, 
    std::vector<CachedIndex_*>& index_ptrs, 
    Index_* counts)
{
    auto ctype = get_class_name(matrix);
    if (ctype != "SVT_SparseMatrix") {
        // Can't be bothered to write a parser for COO_SparseMatrix objects,
        // which are soon-to-be-superceded by SVT_SparseMatrix anyway; so we
        // just forcibly coerce it.
        auto methods_env = Rcpp::Environment::namespace_env("methods");
        Rcpp::Function converter(methods_env["as"]);
        matrix = converter(matrix, Rcpp::CharacterVector::create("SVT_SparseMatrix"));
    }

    bool needs_value = !value_ptrs.empty();
    bool needs_index = !index_ptrs.empty();

    parse_SVT_SparseMatrix(matrix, [&](int c, const auto& curindices, bool all_ones, const auto& curvalues) {
        size_t nnz = curindices.size();

        // Note that non-empty value_ptrs and index_ptrs may be longer than the
        // number of rows/columns in the SVT matrix, due to the reuse of slabs.
        if (row) {
            if (needs_value) {
                if (all_ones) {
                    for (size_t i = 0; i < nnz; ++i) {
                        auto ix = curindices[i];
                        value_ptrs[ix][counts[ix]] = 1;
                    }
                } else {
                    for (size_t i = 0; i < nnz; ++i) {
                        auto ix = curindices[i];
                        value_ptrs[ix][counts[ix]] = curvalues[i];
                    }
                }
            }
            if (needs_index) {
                for (size_t i = 0; i < nnz; ++i) {
                    auto ix = curindices[i];
                    index_ptrs[ix][counts[ix]] = c;
                }
            }
            for (size_t i = 0; i < nnz; ++i) {
                ++(counts[curindices[i]]);
            }

        } else {
            if (needs_value) {
                if (all_ones) {
                    std::fill_n(value_ptrs[c], nnz, 1);
                } else {
                    std::copy(curvalues.begin(), curvalues.end(), value_ptrs[c]);
                }
            }
            if (needs_index) {
                std::copy(curindices.begin(), curindices.end(), index_ptrs[c]);
            }
            counts[c] = nnz;
        }
    });
}
/**
 * @endcond
 */

}

#endif
