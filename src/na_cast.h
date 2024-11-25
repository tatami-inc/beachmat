#ifndef BEACHMAT_NA_CAST_HPP
#define BEACHMAT_NA_CAST_HPP

#include "Rtatami.h"
#include "Rcpp.h"

template<typename Entity_>
bool has_na_integer(const Entity_& mat) {
    for (auto x : mat) {
        if (x == NA_INTEGER) {
            return true;
        }
    }
    return false;
}

template<typename Entity_>
bool has_na_logical(const Entity_& mat) {
    for (auto x : mat) {
        if (x == NA_LOGICAL) {
            return true;
        }
    }
    return false;
}

std::shared_ptr<tatami::NumericMatrix> delayed_cast_na_integer(std::shared_ptr<tatami::NumericMatrix>);

std::shared_ptr<tatami::NumericMatrix> delayed_cast_na_logical(std::shared_ptr<tatami::NumericMatrix>);

#endif
