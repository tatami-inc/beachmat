#include "Rtatami.h"

#include <type_traits>
#include <stdexcept>

#include "na_cast.h"

//[[Rcpp::export(rng=false)]]
SEXP initialize_dense_matrix(Rcpp::RObject raw_x, int nrow, int ncol, bool check_na) {
    auto output = Rtatami::new_BoundNumericMatrix();

    if (raw_x.sexp_type() == INTSXP) {
        Rcpp::IntegerMatrix x(raw_x);
        output->original = x; // Hold reference to avoid GC, in case of allocations to create 'x', e.g., for ALTREP.
        tatami::ArrayView<int> x_view(static_cast<const int*>(x.begin()), x.size());
        output->ptr.reset(new tatami::DenseMatrix<double, int, decltype(x_view)>(nrow, ncol, std::move(x_view), false));
        if (check_na) {
            auto masked = delayed_cast_na_integer(std::move(output->ptr)); 
            output->ptr = std::move(masked);
        }

    } else if (raw_x.sexp_type() == LGLSXP) {
        Rcpp::LogicalMatrix x(raw_x);
        output->original = x;
        tatami::ArrayView<int> x_view(static_cast<const int*>(x.begin()), x.size());
        output->ptr.reset(new tatami::DenseMatrix<double, int, decltype(x_view)>(nrow, ncol, std::move(x_view), false));
        if (check_na) {
            auto masked = delayed_cast_na_logical(std::move(output->ptr)); 
            output->ptr = std::move(masked);
        }

    } else if (raw_x.sexp_type() == REALSXP) {
        Rcpp::NumericMatrix x(raw_x);
        output->original = x; 
        tatami::ArrayView<double> x_view(static_cast<const double*>(x.begin()), x.size());
        output->ptr.reset(new tatami::DenseMatrix<double, int, decltype(x_view)>(nrow, ncol, std::move(x_view), false));

    } else {
        throw std::runtime_error("'x' vector should be integer or real");
    }

    return output;
}

//[[Rcpp::export(rng=false)]]
SEXP initialize_dense_matrix_from_vector(Rcpp::RObject raw_x, int nrow, int ncol, bool check_na) {
    auto output = Rtatami::new_BoundNumericMatrix();

    if (raw_x.sexp_type() == LGLSXP) {
        Rcpp::LogicalVector x(raw_x);
        output->original = x;
        tatami::ArrayView<int> x_view(static_cast<const int*>(x.begin()), x.size());
        output->ptr.reset(new tatami::DenseMatrix<double, int, decltype(x_view)>(nrow, ncol, std::move(x_view), false));
        if (check_na) { 
            auto masked = delayed_cast_na_logical(std::move(output->ptr)); 
            output->ptr = std::move(masked);
        }

    } else if (raw_x.sexp_type() == REALSXP) {
        Rcpp::NumericVector x(raw_x);
        output->original = x; 
        tatami::ArrayView<double> x_view(static_cast<const double*>(x.begin()), x.size());
        output->ptr.reset(new tatami::DenseMatrix<double, int, decltype(x_view)>(nrow, ncol, std::move(x_view), false));

    } else {
        throw std::runtime_error("'x' vector should be integer or real");
    }

    return output;
}
