#include "tatamize.h"
#include "tatami/tatami.hpp"
#include "Rcpp.h"

#include <cstdint>
#include <limits>
#include <type_traits>

template<class V>
struct RcppVectorPlus {
    RcppVectorPlus(V in) : vec(std::move(in)) {}
    typedef typename std::remove_reference<decltype(V()[0])>::type T;

    const T* data() const {
        return static_cast<const T*>(vec.begin());
    }

    auto begin() const {
        return vec.begin();
    }

    auto end() const {
        return vec.end();
    }

    size_t size() const {
        return vec.size();
    }

    T operator[](size_t i) const {
        return vec[i];
    }
private:
    V vec;
};

template<class XVector, class IVector, class PVector>
SEXP create_matrix_nocopy(XVector x, IVector i, PVector p, int nrow, int ncol, bool byrow) {
    RcppVectorPlus<XVector> x_(x);
    RcppVectorPlus<IVector> i_(i);

    if (byrow) {
        typedef tatami::CompressedSparseMatrix<true, double, int, decltype(x_), decltype(i_), decltype(p)> SparseMat;
        return new_MatrixChan(new SparseMat(nrow, ncol, std::move(x_), std::move(i_), std::move(p), false));
    } else {
        typedef tatami::CompressedSparseMatrix<false, double, int, decltype(x_), decltype(i_), decltype(p)> SparseMat;
        return new_MatrixChan(new SparseMat(nrow, ncol, std::move(x_), std::move(i_), std::move(p), false));
    }
}

//[[Rcpp::export(rng=false)]]
SEXP initialize_sparse_matrix(Rcpp::RObject x, Rcpp::RObject i, Rcpp::RObject p, int nrow, int ncol, bool byrow) {
    if (p.sexp_type() != INTSXP) {
        throw std::runtime_error("'p' vector should be integer");
    }
    Rcpp::IntegerVector p_(p);

    if (i.sexp_type() != INTSXP) {
        throw std::runtime_error("'i' vector should be integer");
    }
    Rcpp::IntegerVector i_(i);

    if (x.sexp_type() == INTSXP) {
        Rcpp::IntegerVector x_(x);
        return create_matrix_nocopy(std::move(x_), std::move(i_), std::move(p_), nrow, ncol, byrow);
    } else if (x.sexp_type() != REALSXP) {
        throw std::runtime_error("'x' vector should be integer or real");
    }

    Rcpp::NumericVector x_(x);
    return create_matrix_nocopy(std::move(x_), std::move(i_), std::move(p_), nrow, ncol, byrow);
}
