#ifndef UTILS_H
#define UTILS_H

#include "Rcpp.h"

#include <type_traits>

template<typename Input_>
using I = std::remove_cv_t<std::remove_reference_t<Input_> >;

template<typename Row_, typename Column_>
Rcpp::NumericMatrix create_output_matrix(Row_ nrow, Column_ ncol) {
    return Rcpp::NumericMatrix(
        I<decltype(std::declval<Rcpp::NumericMatrix>().rows())>(nrow),
        I<decltype(std::declval<Rcpp::NumericMatrix>().cols())>(ncol)
    );
}

#endif
