#ifndef TATAMI_R_SIMPLEMATRIX_HPP
#define TATAMI_R_SIMPLEMATRIX_HPP

#include "utils.hpp"
#include "tatami/tatami.hpp"
#include "tatami/utils/ArrayView.hpp"

namespace tatami_r { 

template<typename Data_ = double, typename Index_ = int, class InputObject_>
Parsed<Data_, Index_> parse_simple_matrix_internal(const InputObject_& y) {
    Parsed<Data_, Index_> output;

    typedef typename std::remove_const<typename std::remove_reference<decltype(y[0])>::type>::type Value_;
    tatami::ArrayView view(static_cast<const Value_*>(y.begin()), y.size());
    output.matrix.reset(new tatami::DenseColumnMatrix<double, int, decltype(view)>(y.rows(), y.cols(), std::move(view)));

    output.contents = Rcpp::List::create(y);
    return output;
}

template<typename Data_ = double, typename Index_ = int>
Parsed<Data_, Index_> parse_simple_matrix(const Rcpp::RObject& seed) {
    Parsed<Data_, Index_> output;

    if (seed.sexp_type() == REALSXP) {
        Rcpp::NumericMatrix y(seed);
        output = parse_simple_matrix_internal<Data_, Index_>(y);
    } else if (seed.sexp_type() == INTSXP) {
        Rcpp::IntegerMatrix y(seed);
        output = parse_simple_matrix_internal<Data_, Index_>(y);
    } else if (seed.sexp_type() == LGLSXP) {
        Rcpp::LogicalMatrix y(seed);
        output = parse_simple_matrix_internal<Data_, Index_>(y);
    }

    return output;
}

}

#endif
