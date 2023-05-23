#include "Rtatami.h"

#include <cstdint>
#include <limits>
#include <type_traits>

//[[Rcpp::export(rng=false)]]
SEXP initialize_dense_matrix(Rcpp::RObject raw_x, int nrow, int ncol) {
    auto output = Rtatami::new_BoundNumericMatrix();
    output->original = raw_x; // Hold reference to avoid GC.

    if (raw_x.sexp_type() == INTSXP) {
        Rcpp::IntegerVector x(raw_x);
        tatami::ArrayView<int> x_view(static_cast<const int*>(x.begin()), x.size());
        output->ptr.reset(new tatami::DenseColumnMatrix<double, int, decltype(x_view)>(nrow, ncol, std::move(x_view)));
    } else if (raw_x.sexp_type() == LGLSXP) {
        Rcpp::LogicalVector x(raw_x);
        tatami::ArrayView<int> x_view(static_cast<const int*>(x.begin()), x.size());
        output->ptr.reset(new tatami::DenseColumnMatrix<double, int, decltype(x_view)>(nrow, ncol, std::move(x_view)));
    } else if (raw_x.sexp_type() == REALSXP) {
        Rcpp::NumericVector x(raw_x);
        tatami::ArrayView<double> x_view(static_cast<const double*>(x.begin()), x.size());
        output->ptr.reset(new tatami::DenseColumnMatrix<double, int, decltype(x_view)>(nrow, ncol, std::move(x_view)));
    } else {
        throw std::runtime_error("'x' vector should be integer or real");
    }

    return output;
}
