#include "Rtatami.h"

#include <cstdint>
#include <limits>
#include <type_traits>

template<typename T_, typename XVector_>
tatami::NumericMatrix* store_sparse_matrix(XVector_ x, Rcpp::IntegerVector i, Rcpp::IntegerVector p, int nrow, int ncol, bool byrow) {
    tatami::ArrayView<T_> x_view(static_cast<const T_*>(x.begin()), x.size());
    tatami::ArrayView<int> i_view(static_cast<const int*>(i.begin()), i.size());
    tatami::ArrayView<int> p_view(static_cast<const int*>(p.begin()), p.size());

    if (byrow) {
        typedef tatami::CompressedSparseMatrix<true, double, int, decltype(x_view), decltype(i_view), decltype(p_view)> SparseMat;
        return new SparseMat(nrow, ncol, std::move(x_view), std::move(i_view), std::move(p_view), false);
    } else {
        typedef tatami::CompressedSparseMatrix<false, double, int, decltype(x_view), decltype(i_view), decltype(p_view)> SparseMat;
        return new SparseMat(nrow, ncol, std::move(x_view), std::move(i_view), std::move(p_view), false);
    }
}

//[[Rcpp::export(rng=false)]]
SEXP initialize_sparse_matrix(Rcpp::RObject raw_x, Rcpp::RObject raw_i, Rcpp::RObject raw_p, int nrow, int ncol, bool byrow) {
    auto output = Rtatami::new_BoundNumericMatrix();
    output->original = Rcpp::List::create(raw_x, raw_i, raw_p); // holding references to all R objects used here.

    if (raw_p.sexp_type() != INTSXP) {
        throw std::runtime_error("'p' vector should be integer");
    }
    Rcpp::IntegerVector p(raw_p);

    if (raw_i.sexp_type() != INTSXP) {
        throw std::runtime_error("'i' vector should be integer");
    }
    Rcpp::IntegerVector i(raw_i);

    if (raw_x.sexp_type() == LGLSXP) {
        Rcpp::LogicalVector x(raw_x);
        output->ptr.reset(store_sparse_matrix<int>(std::move(x), std::move(i), std::move(p), nrow, ncol, byrow));
    } else if (raw_x.sexp_type() == REALSXP) {
        Rcpp::NumericVector x(raw_x);
        output->ptr.reset(store_sparse_matrix<double>(std::move(x), std::move(i), std::move(p), nrow, ncol, byrow));
    } else {
        throw std::runtime_error("'x' vector should be integer or real");
    }

    return output;
}
