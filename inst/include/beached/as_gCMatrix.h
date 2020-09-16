#ifndef BEACHMAT_AS_GCMATRIX_H
#define BEACHMAT_AS_GCMATRIX_H

#include "Rcpp.h"
#include <map>

namespace beachmat {

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

template <class V>
inline Rcpp::RObject as_gCMatrix (int nr, int nc, V x, Rcpp::IntegerVector i, Rcpp::IntegerVector p) {
    auto mat = generate_gCMatrix<V>();
    mat.slot("Dim") = Rcpp::IntegerVector::create(nr, nc);
    mat.slot("p")=p;
    mat.slot("i")=i;
    mat.slot("x")=x;
    return SEXP(mat);
}

}

#endif
