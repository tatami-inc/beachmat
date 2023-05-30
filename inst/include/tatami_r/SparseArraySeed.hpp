#ifndef TATAMI_R_SPARSEARRAYSEED_HPP
#define TATAMI_R_SPARSEARRAYSEED_HPP

#include "utils.hpp"
#include "tatami/tatami.hpp"
#include <type_traits>

namespace tatami_r { 

template<typename Data_ = double, typename Index_ = int, class InputObject_>
Parsed<Data_, Index_> parse_SparseArraySeed_internal(Rcpp::RObject seed, InputObject_ val, bool prefer_csr) {
    auto dims = parse_dims(seed.slot("dim"));
    int NR = dims.first;
    int NC = dims.second;

    Rcpp::IntegerMatrix temp_i(Rcpp::RObject(seed.slot("nzindex")));
    if (temp_i.ncol() != 2) {
        auto ctype = get_class_name(seed);
        throw std::runtime_error(std::string("'nzindex' slot in a ") + ctype + " object should have two columns"); 
    }

    const size_t nnz = temp_i.nrow();
    if (nnz != val.size()) {
        auto ctype = get_class_name(seed);
        throw std::runtime_error(std::string("incompatible 'nzindex' and 'nzdata' lengths in a ") + ctype + " object"); 
    }

    auto row_indices = temp_i.column(0);
    auto col_indices = temp_i.column(1);

    // Checking if it's already sorted, which allows us to skip our own re-sorting.
    bool okay_csc = true;
    bool okay_csr = true;
    if (nnz) {
        auto rowIt = row_indices.begin();
        auto colIt = col_indices.begin();

        auto check_index = [&](int r, int c) -> void {
            if (r <= 0 || r > NR || c <= 0 || c > NC) {
                auto ctype = get_class_name(seed);
                throw std::runtime_error(std::string("'nzindex' out of bounds in a ") + ctype + " object");
            }
        };

        auto lastR = *rowIt;
        auto lastC = *colIt;
        check_index(lastR, lastC);

        for (size_t v = 1; v < nnz; ++v) {
            auto nextR = *(++rowIt);
            auto nextC = *(++colIt);
            check_index(nextR, nextC);

            if (okay_csc) {
                if (lastC > nextC || (lastC == nextC && lastR > nextR)) {
                    okay_csc = false;
                }
            }
            if (okay_csr) {
                if (lastR > nextR || (lastR == nextR && lastC > nextC)) {
                    okay_csr = false;
                }
            }

            lastC = nextC;
            lastR = nextR;
        }
    }

    typedef typename std::remove_const<typename std::remove_reference<decltype(val[0])>::type>::type Value_;

    Parsed<Data_, Index_> output;
    if (okay_csc) {
        std::vector<int> i(row_indices.begin(), row_indices.end());
        for (auto& x : i) {
            --x;
        }

        std::vector<size_t> p(NC + 1);
        auto colIt = col_indices.begin();
        auto colStart = colIt;
        for (int c = 1; c <= NC; ++c) {
            // Technically it should be *colIt <= c+1 to get to 1-based
            // indices, but this cancels out with a -1 because we want
            // everything up to the _last_ column.
            while (colIt != col_indices.end() && *colIt <= c) { 
                ++colIt;
            }
            p[c] = colIt - colStart;
        }

        tatami::ArrayView<Value_> vview(static_cast<const Value_*>(val.begin()), val.size());
        output.contents = std::move(val);

        output.matrix.reset(
            new tatami::CompressedSparseMatrix<false, Data_, Index_, decltype(vview), decltype(i), decltype(p)>(
                NR, 
                NC, 
                std::move(vview), 
                std::move(i), 
                std::move(p), 
                false
            )
        );

    } else if (okay_csr) {
        std::vector<int> j(col_indices.begin(), col_indices.end());
        for (auto& x : j) {
            --x;
        }

        std::vector<size_t> p(NR + 1);
        auto rowIt = row_indices.begin();
        auto rowStart = rowIt;
        for (int r = 1; r <= NR; ++r) {
            // Technically it should be *rowIt <= r+1 to get to 1-based
            // indices, but this cancels out with a -1 because we want
            // everything up to the _last_ rowumn.
            while (rowIt != row_indices.end() && *rowIt <= r) { 
                ++rowIt;
            }
            p[r] = rowIt - rowStart;
        }

        tatami::ArrayView<Value_> vview(static_cast<const Value_*>(val.begin()), val.size());
        output.contents = std::move(val);

        output.matrix.reset(
            new tatami::CompressedSparseMatrix<true, Data_, Index_, decltype(vview), decltype(j), decltype(p)>(
                NR, 
                NC, 
                std::move(vview), 
                std::move(j), 
                std::move(p), 
                false
            )
        );

    } else {
        std::vector<Value_> v(val.begin(), val.end());
        std::vector<int> i(row_indices.begin(), row_indices.end());
        for (auto& x : i) {
            --x;
        }
        std::vector<int> j(col_indices.begin(), col_indices.end());
        for (auto& x : j) {
            --x;
        }
        output.contents = R_NilValue; // no need for anything from the original object.

        // If we need to sort anyway, we might as well prefer to the layout
        // that's more amenable to the UnknownMatrix's extraction pattern.
        if (prefer_csr) {
            auto p = tatami::compress_sparse_triplets<true>(NR, NC, v, i, j);
            output.matrix.reset(
                new tatami::CompressedSparseMatrix<true, Data_, Index_, decltype(v), decltype(j), decltype(p)>(
                    NR, 
                    NC, 
                    std::move(v), 
                    std::move(j), 
                    std::move(p), 
                    false 
                )
            );

        } else {
            auto p = tatami::compress_sparse_triplets<false>(NR, NC, v, i, j);
            output.matrix.reset(
                new tatami::CompressedSparseMatrix<false, Data_, Index_, decltype(v), decltype(i), decltype(p)>(
                    NR, 
                    NC, 
                    std::move(v), 
                    std::move(i), 
                    std::move(p), 
                    false
                )
            );
        }
    }

    return output;
}

template<typename Data_ = double, typename Index_ = int>
Parsed<Data_, Index_> parse_SparseArraySeed(Rcpp::RObject seed, bool prefer_csr) {
    Rcpp::RObject vals(seed.slot("nzdata"));

    Parsed<Data_, Index_> output;
    if (vals.sexp_type() == REALSXP) {
        output = parse_SparseArraySeed_internal<Data_, Index_>(seed, Rcpp::NumericVector(vals), prefer_csr);
    } else if (vals.sexp_type() == INTSXP) {
        output = parse_SparseArraySeed_internal<Data_, Index_>(seed, Rcpp::IntegerVector(vals), prefer_csr);
    } else if (vals.sexp_type() == LGLSXP) {
        output = parse_SparseArraySeed_internal<Data_, Index_>(seed, Rcpp::LogicalVector(vals), prefer_csr);
    }

    return output;
}

}

#endif
