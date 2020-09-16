#ifndef BEACHMAT_AS_GCMATRIX_H
#define BEACHMAT_AS_GCMATRIX_H

/**
 * @file as_gCMatrix.h
 *
 * Functions to create `*gCMatrix` instances.
 */

#include "Rcpp.h"
#include <map>

namespace beachmat {

/** 
 * @internal
 *
 * Generate a `*gCMatrix` of the appropriate class for the `Rcpp::Vector` class.
 *
 * @tparam V An `Rcpp::Vector` class.
 * Currently only `Rcpp::NumericVector` and `Rcpp::LogicalVector` are supported.
 *
 * @return An `Rcpp::S4` object of the appropriate `*gCMatrix` class for `V`,
 * either `dgCMatrix` or `lgCMatrix` for double-precision and logical data, respectively.
 * All other `V` will trigger a compile-time error.
 */
template <class V>
inline Rcpp::S4 generate_gCMatrix () { // could also use = delete here.
    static_assert(sizeof(V)==0, "unsupported specialization of generate_gCMatrix");
}

template <>
inline Rcpp::S4 generate_gCMatrix<Rcpp::LogicalVector> () {
    return Rcpp::S4("lgCMatrix");
}

template <>
inline Rcpp::S4 generate_gCMatrix<Rcpp::NumericVector> () {
    return Rcpp::S4("dgCMatrix");
}

/**
 * Create a `*gCMatrix` from a triplet map of non-zero entries.
 * Best used when the number of non-zero entries is not known in advance.
 *
 * @tparam V An `Rcpp::Vector` class.
 * Only `Rcpp::NumericVector` and `Rcpp::LogicalVector` are supported.
 * @tparam T Type of data to be stored in the `V` class.
 *
 * @param nr Number of rows.
 * @param nc Number of columns.
 * @param holder Triplet-formatted information for all non-zero entries.
 * The key should contain the zero-based column index (first) and the zero-based row index (second).
 * The value should contain the non-zero value itself.
 *
 * @return A `dgCMatrix` or `lgCMatrix` instance (depending on `V`) containing all entries in `holder`.
 */
template <class V, typename T = typename V::stored_type>
inline Rcpp::RObject as_gCMatrix (int nr, int nc, const std::map<std::pair<int, int>, T>& holder) {
    auto mat = generate_gCMatrix<V>();
    mat.slot("Dim") = Rcpp::IntegerVector::create(nr, nc);

    // Setting the 'p' vector.
    Rcpp::IntegerVector p(nc + 1, 0);

    auto hIt = holder.begin();
    for (int c = 0; c <= nc; ++c) {
        while (hIt != holder.end() && (hIt->first).first <= c) {
            ++hIt;
        }
        p[c] = (hIt - holder.begin());
    }
    mat.slot("p")=p;

    // Setting 'i' and 'x'.
    size_t total_size = holder.size();
    Rcpp::IntegerVector i(total_size);
    V x(total_size);

    auto xIt=x.begin();
    auto iIt=i.begin();
    for (auto hIt = holder.begin(); hIt != holder.end(); ++hIt, ++xIt, ++iIt) {
        (*xIt) = (hIt->second);
        (*iIt) = (hIt->first).second;
    }

    mat.slot("i")=i;
    mat.slot("x")=x;

    return SEXP(mat);
}

/**
 * Create a `*gCMatrix` containing modified values from an existing `*gCMatrix` in the same order.
 *
 * @tparam V An `Rcpp::Vector` class.
 * Only `Rcpp::NumericVector` and `Rcpp::LogicalVector` are supported.
 * @tparam T Type of data to be stored in the `V` class.
 *
 * @param old A `*gCMatrix` object.
 * This need not be of a type corresponding to `V`.
 * @param x A vector of non-zero values to use to replace the existing `x` slot in `old`.
 * This should correspond to the same positions and ordering of the non-zero values already in `old`.
 *
 * @return A `dgCMatrix` or `lgCMatrix` instance (depending on `V`) 
 * containing indexing information from `old` with the values in `x`.
 */
template <class V>
inline Rcpp::RObject as_gCMatrix (Rcpp::RObject old, V x) {
    auto mat = generate_gCMatrix<V>();

    mat.slot("Dim") = get_safe_slot(old, "Dim");
    mat.slot("p") = get_safe_slot(old, "p");
    mat.slot("i") = get_safe_slot(old, "i");
    mat.slot("x") = x;
    return SEXP(mat);
}

}

#endif
