#include "Rtatami.h"
#include "Rcpp.h"
#include "tatami/tatami.hpp"
#include "Rmath.h"

#include <memory>
#include <string>
#include <stdexcept>

//[[Rcpp::export(rng=false)]]
SEXP apply_delayed_log(SEXP raw_input, double base) {
    Rtatami::BoundNumericPointer input(raw_input);
    auto output = Rtatami::new_BoundNumericMatrix();
    output->ptr.reset(new tatami::DelayedUnaryIsometricOperation<double, double, int>(input->ptr, std::make_shared<tatami::DelayedUnaryIsometricLogHelper<double, double, int, double> >(base)));
    output->original = input->original; // copying the reference to propagate GC protection.
    return output;
}

//[[Rcpp::export(rng=false)]]
SEXP apply_delayed_unary_math(SEXP raw_input, const std::string& op) {
    Rtatami::BoundNumericPointer input(raw_input);
    auto iptr = input->ptr;

    std::shared_ptr<tatami::DelayedUnaryIsometricOperationHelper<double, double, int> > opptr;
    if (op == "abs") {
        opptr.reset(new tatami::DelayedUnaryIsometricAbsHelper<double, double, int>);
    } else if (op == "sign") {
        opptr.reset(new tatami::DelayedUnaryIsometricSignHelper<double, double, int>);
    } else if (op == "sqrt") {
        opptr.reset(new tatami::DelayedUnaryIsometricSqrtHelper<double, double, int>);
    } else if (op == "floor") {
        opptr.reset(new tatami::DelayedUnaryIsometricFloorHelper<double, double, int>);
    } else if (op == "ceiling") {
        opptr.reset(new tatami::DelayedUnaryIsometricCeilingHelper<double, double, int>);
    } else if (op == "trunc") {
        opptr.reset(new tatami::DelayedUnaryIsometricTruncHelper<double, double, int>);
    } else if (op == "exp") {
        opptr.reset(new tatami::DelayedUnaryIsometricExpHelper<double, double, int>);
    } else if (op == "expm1") {
        opptr.reset(new tatami::DelayedUnaryIsometricExpm1Helper<double, double, int>);
    } else if (op == "log1p") {
        opptr.reset(new tatami::DelayedUnaryIsometricLog1pHelper<double, double, int, double>);
    } else if (op == "cos") {
        opptr.reset(new tatami::DelayedUnaryIsometricCosHelper<double, double, int>);
    } else if (op == "sin") {
        opptr.reset(new tatami::DelayedUnaryIsometricSinHelper<double, double, int>);
    } else if (op == "tan") {
        opptr.reset(new tatami::DelayedUnaryIsometricTanHelper<double, double, int>);
    } else if (op == "cospi" || op == "sinpi" || op == "tanpi") {
        auto tmp_ptr = std::make_shared<tatami::DelayedUnaryIsometricOperation<double, double, int> >(
            std::move(iptr),
            std::make_shared<tatami::DelayedUnaryIsometricMultiplyScalarHelper<double, double, int, double> >(M_PI)
        );
        iptr = std::move(tmp_ptr);
        if (op == "cospi") {
            opptr.reset(new tatami::DelayedUnaryIsometricCosHelper<double, double, int>);
        } else if (op == "sinpi") {
            opptr.reset(new tatami::DelayedUnaryIsometricSinHelper<double, double, int>);
        } else {
            opptr.reset(new tatami::DelayedUnaryIsometricTanHelper<double, double, int>);
        }
    } else if (op == "acos") {
        opptr.reset(new tatami::DelayedUnaryIsometricAcosHelper<double, double, int>());
    } else if (op == "asin") {
        opptr.reset(new tatami::DelayedUnaryIsometricAsinHelper<double, double, int>());
    } else if (op == "atan") {
        opptr.reset(new tatami::DelayedUnaryIsometricAtanHelper<double, double, int>());
    } else if (op == "cosh") {
        opptr.reset(new tatami::DelayedUnaryIsometricCoshHelper<double, double, int>());
    } else if (op == "sinh") {
        opptr.reset(new tatami::DelayedUnaryIsometricSinhHelper<double, double, int>());
    } else if (op == "tanh") {
        opptr.reset(new tatami::DelayedUnaryIsometricTanhHelper<double, double, int>());
    } else if (op == "acosh") {
        opptr.reset(new tatami::DelayedUnaryIsometricAcoshHelper<double, double, int>());
    } else if (op == "asinh") {
        opptr.reset(new tatami::DelayedUnaryIsometricAsinhHelper<double, double, int>());
    } else if (op == "atanh") {
        opptr.reset(new tatami::DelayedUnaryIsometricAtanhHelper<double, double, int>());
    } else if (op == "lgamma") {
        opptr.reset(new tatami::DelayedUnaryIsometricLgammaHelper<double, double, int>());
    } else if (op == "gamma") {
        opptr.reset(new tatami::DelayedUnaryIsometricGammaHelper<double, double, int>());
    } else {
        return(R_NilValue);
    }

    auto output = Rtatami::new_BoundNumericMatrix();
    output->ptr.reset(new tatami::DelayedUnaryIsometricOperation<double, double, int>(std::move(iptr), std::move(opptr)));
    output->original = input->original; // copying the reference to propagate GC protection.
    return output;
}

//[[Rcpp::export(rng=false)]]
SEXP apply_delayed_round(SEXP raw_input) {
    Rtatami::BoundNumericPointer input(raw_input);
    auto output = Rtatami::new_BoundNumericMatrix();
    output->ptr.reset(new tatami::DelayedUnaryIsometricOperation<double, double, int>(input->ptr, std::make_shared<tatami::DelayedUnaryIsometricRoundHelper<double, double, int> >()));
    output->original = input->original; // copying the reference to propagate GC protection.
    return output;
}
