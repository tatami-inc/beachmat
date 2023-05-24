#ifndef TATAMI_R_SPARSEARRAYSEED_HPP
#define TATAMI_R_SPARSEARRAYSEED_HPP

#include "utils.hpp"
#include "tatami/tatami.hpp"
#include <type_traits>

namespace tatami_r { 

template<typename Data_ = double, typename Index_ = int, class InputObject_>
Parsed<Data_, Index_> parse_SparseArraySeed_internal(Rcpp::RObject seed, const InputObject_& val) {
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

    std::vector<int> i(nnz);
    std::vector<size_t> p(NC + 1);
    typedef typename std::remove_const<typename std::remove_reference<decltype(val[0])>::type>::type Value_;
    std::vector<Value_> v(val.begin(), val.end());

    if (nnz) {
        auto row_indices=temp_i.column(0);
        auto col_indices=temp_i.column(1);

        auto iIt = i.begin();
        for (auto subi : row_indices) { 
            *iIt = subi - 1;
            ++iIt;
        }

        // Checking if it's already sorted.
        bool okay=true;
        {
            auto rowIt=row_indices.begin();
            auto colIt=col_indices.begin();

            for (size_t v = 0; v < nnz; ++v) {
                auto lastR=*rowIt;
                auto lastC=*colIt;
                if (lastR <= 0 || lastR > NR || lastC <= 0 || lastC > NC) {
                    auto ctype = get_class_name(seed);
                    throw std::runtime_error(std::string("'nzindex' out of bounds in a ") + ctype + " object");
                }

                if (okay && v < nnz - 1) {
                    auto nextR=*(++rowIt);
                    auto nextC=*(++colIt);
                    if (lastC > nextC || (lastC == nextC && lastR > nextR)) {
                        okay = false;
                    }
                }
            }
        }

        if (okay) {
            p.resize(NC + 1);
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

        } else {
            std::vector<int> j(nnz);
            auto jIt = j.begin();
            for (auto subj : row_indices) { 
                *jIt = subj - 1;
                ++jIt;
            }
            p = tatami::compress_sparse_triplets<false>(NR, NC, v, i, j);
        }
    }

    Parsed<Data_, Index_> output;
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
    return output;
}

template<typename Data_ = double, typename Index_ = int>
Parsed<Data_, Index_> parse_SparseArraySeed(Rcpp::RObject seed) {
    Rcpp::RObject vals(seed.slot("nzdata"));

    Parsed<Data_, Index_> output;
    if (vals.sexp_type() == REALSXP) {
        Rcpp::NumericVector y(vals);
        output = parse_SparseArraySeed_internal<Data_, Index_>(seed, y);
    } else if (vals.sexp_type() == INTSXP) {
        Rcpp::IntegerVector y(vals);
        output = parse_SparseArraySeed_internal<Data_, Index_>(seed, y);
    } else if (vals.sexp_type() == LGLSXP) {
        Rcpp::LogicalVector y(vals);
        output = parse_SparseArraySeed_internal<Data_, Index_>(seed, y);
    }

    return output;
}

}

#endif
