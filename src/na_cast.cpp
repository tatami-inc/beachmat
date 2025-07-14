#include "na_cast.h"

#include <memory>

std::shared_ptr<tatami::NumericMatrix> delayed_cast_na_integer(std::shared_ptr<tatami::NumericMatrix> seed) {
    return std::make_shared<tatami::DelayedUnaryIsometricOperation<double, double, int> >(
        std::move(seed),
        std::make_shared<tatami::DelayedUnaryIsometricSubstituteScalarHelper<tatami::CompareOperation::EQUAL, double, double, int, double> >(NA_INTEGER, NA_REAL)
    );
}

std::shared_ptr<tatami::NumericMatrix> delayed_cast_na_logical(std::shared_ptr<tatami::NumericMatrix> seed) {
    return std::make_shared<tatami::DelayedUnaryIsometricOperation<double, double, int> >(
        std::move(seed),
        std::make_shared<tatami::DelayedUnaryIsometricSubstituteScalarHelper<tatami::CompareOperation::EQUAL, double, double, int, double> >(NA_LOGICAL, NA_REAL)
    );
}
