#ifndef BEACHMAT_NA_CAST_HPP
#define BEACHMAT_NA_CAST_HPP

#include "Rtatami.h"
#include "Rcpp.h"

std::shared_ptr<tatami::NumericMatrix> delayed_cast_na_integer(std::shared_ptr<tatami::NumericMatrix>);

std::shared_ptr<tatami::NumericMatrix> delayed_cast_na_logical(std::shared_ptr<tatami::NumericMatrix>);

#endif
