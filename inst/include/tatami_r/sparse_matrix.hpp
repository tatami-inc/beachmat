#ifndef TATAMI_R_SPARSE_MATRIX_HPP
#define TATAMI_R_SPARSE_MATRIX_HPP

#include "utils.hpp"
#include "tatami/tatami.hpp"
#include <type_traits>

namespace tatami_r { 

template<class InputObject_, SEXPTYPE desired_sexp_, typename InputValue_, typename CachedValue_, typename CachedIndex_, typename Index_>
void parse_sparse_matrix_internal(
    Rcpp::RObject seed, 
    bool row,
    std::vector<CachedValue_*>& value_ptrs, 
    std::vector<CachedIndex_*>& index_ptrs, 
    Index_* counts)
{
    Rcpp::RObject raw_svt = seed.slot("SVT");
    if (raw_svt == R_NilValue) {
        return;
    }

    Rcpp::IntegerVector svt_version(seed.slot(".svt_version"));
    if (svt_version.size() != 1) {
        throw std::runtime_error("'.svt_version' should be an integer scalar");
    }
    int version = svt_version[0];
    int index_x = (version == 0 ? 0 : 1);
    int value_x = (version == 0 ? 1 : 0);

    Rcpp::List svt(raw_svt);
    int NC = svt.size();
    bool needs_value = !value_ptrs.empty();
    bool needs_index = !index_ptrs.empty();

    // Note that non-empty value_ptrs and index_ptrs may be longer than the
    // number of rows/columns in the SVT matrix, due to the reuse of slabs.

    for (int c = 0; c < NC; ++c) {
        Rcpp::RObject raw_inner(svt[c]);
        if (raw_inner == R_NilValue) {
            continue;
        }

        Rcpp::List inner(raw_inner);
        if (inner.size() != 2) {
            auto ctype = get_class_name(seed);
            throw std::runtime_error("each entry of the 'SVT' slot of a " + ctype + " object should be a list of length 2 or NULL");
        }

        // Verify type to ensure that we're not making a view on a temporary array.
        Rcpp::RObject raw_indices = inner[index_x];
        if (raw_indices.sexp_type() != INTSXP) {
            auto ctype = get_class_name(seed);
            throw std::runtime_error("indices of each element of the 'SVT' slot in a " + ctype + " object should be an integer vector");
        }
        Rcpp::IntegerVector curindices(raw_indices);

        InputObject_ curvalues;
        Rcpp::RObject raw_values(inner[value_x]);
        if (raw_values != R_NilValue) {
            if (raw_values.sexp_type() != desired_sexp_) {
                auto ctype = get_class_name(seed);
                throw std::runtime_error("value vector of an element of the 'SVT' slot in a " + ctype + " object has an unexpected type");
            }
            curvalues = InputObject_(raw_values);
        } else {
            curvalues = InputObject_(curindices.size());
            std::fill(curvalues.begin(), curvalues.end(), 1);
        }

        size_t nnz = curvalues.size();
        if (nnz != static_cast<size_t>(curindices.size())) {
            auto ctype = get_class_name(seed);
            throw std::runtime_error("both vectors of an element of the 'SVT' slot in a " + ctype + " object should have the same length");
        }

        if (row) {
            if (needs_value) {
                for (size_t i = 0; i < nnz; ++i) {
                    auto ix = curindices[i];
                    value_ptrs[ix][counts[ix]] = curvalues[i];
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
                std::copy(curvalues.begin(), curvalues.end(), value_ptrs[c]);
            }
            if (needs_index) {
                std::copy(curindices.begin(), curindices.end(), index_ptrs[c]);
            }
            counts[c] = nnz;
        }
    }
}

template<typename CachedValue_, typename CachedIndex_, typename Index_>
void parse_sparse_matrix(
    Rcpp::RObject seed,
    bool row,
    std::vector<CachedValue_*>& value_ptrs, 
    std::vector<CachedIndex_*>& index_ptrs, 
    Index_* counts)
{
    auto ctype = get_class_name(seed);
    if (ctype != "SVT_SparseMatrix") {
        // Can't be bothered to write a parser for COO_SparseMatrix objects,
        // which are soon-to-be-superceded by SVT_SparseMatrix anyway; so we
        // just forcibly coerce it.
        auto methods_env = Rcpp::Environment::namespace_env("methods");
        Rcpp::Function converter(methods_env["as"]);
        seed = converter(seed, Rcpp::CharacterVector::create("SVT_SparseMatrix"));
    }

    std::string type = Rcpp::as<std::string>(seed.slot("type"));
    if (type == "double") {
        parse_sparse_matrix_internal<Rcpp::NumericVector, REALSXP, double>(seed, row, value_ptrs, index_ptrs, counts);
    } else if (type == "integer") {
        parse_sparse_matrix_internal<Rcpp::IntegerVector, INTSXP, int>(seed, row, value_ptrs, index_ptrs, counts);
    } else if (type == "logical") {
        parse_sparse_matrix_internal<Rcpp::LogicalVector, LGLSXP, int>(seed, row, value_ptrs, index_ptrs, counts);
    } else {
        throw std::runtime_error("unsupported type '" + type + "' for a " + ctype);
    }
}

}

#endif
