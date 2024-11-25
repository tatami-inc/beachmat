#include "na_cast.h"

#include <memory>

std::shared_ptr<tatami::NumericMatrix> delayed_cast_na_integer(std::shared_ptr<tatami::NumericMatrix> seed) {
    return tatami::make_DelayedUnaryIsometricOperation(std::move(seed), tatami::DelayedUnaryIsometricSubstituteScalar<tatami::CompareOperation::EQUAL, double>(NA_INTEGER, NA_REAL));
}

std::shared_ptr<tatami::NumericMatrix> delayed_cast_na_logical(std::shared_ptr<tatami::NumericMatrix> seed) {
    return tatami::make_DelayedUnaryIsometricOperation(std::move(seed), tatami::DelayedUnaryIsometricSubstituteScalar<tatami::CompareOperation::EQUAL, double>(NA_LOGICAL, NA_REAL));
}
