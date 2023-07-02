#include "Rtatami.h"
#include "Rcpp.h"
#include "tatami/tatami.hpp"
#include "Rmath.h"

//[[Rcpp::export(rng=false)]]
SEXP apply_delayed_log(SEXP raw_input, double base) {
    Rtatami::BoundNumericPointer input(raw_input);
    auto output = Rtatami::new_BoundNumericMatrix();
    output->ptr = tatami::make_DelayedUnaryIsometricOp(input->ptr, tatami::DelayedLogHelper(base));
    output->original = input->original; // copying the reference to propagate GC protection.
    return output;
}

//[[Rcpp::export(rng=false)]]
SEXP apply_delayed_unary_math(SEXP raw_input, const std::string& op) {
    Rtatami::BoundNumericPointer input(raw_input);
    auto output = Rtatami::new_BoundNumericMatrix();

    if (op == "abs") {
        output->ptr = tatami::make_DelayedUnaryIsometricOp(input->ptr, tatami::DelayedAbsHelper<>());
    } else if (op == "sign") {
        output->ptr = tatami::make_DelayedUnaryIsometricOp(input->ptr, tatami::DelayedSignHelper<>());
    } else if (op == "sqrt") {
        output->ptr = tatami::make_DelayedUnaryIsometricOp(input->ptr, tatami::DelayedSqrtHelper<>());
    } else if (op == "floor") {
        output->ptr = tatami::make_DelayedUnaryIsometricOp(input->ptr, tatami::DelayedFloorHelper<>());
    } else if (op == "ceiling") {
        output->ptr = tatami::make_DelayedUnaryIsometricOp(input->ptr, tatami::DelayedCeilingHelper<>());
    } else if (op == "trunc") {
        output->ptr = tatami::make_DelayedUnaryIsometricOp(input->ptr, tatami::DelayedTruncHelper<>());

    } else if (op == "exp") {
        output->ptr = tatami::make_DelayedUnaryIsometricOp(input->ptr, tatami::DelayedExpHelper<>());
    } else if (op == "expm1") {
        output->ptr = tatami::make_DelayedUnaryIsometricOp(input->ptr, tatami::DelayedExpm1Helper<>());
    } else if (op == "log1p") {
        output->ptr = tatami::make_DelayedUnaryIsometricOp(input->ptr, tatami::DelayedLog1pHelper<>());
    } else if (op == "cos") {
        output->ptr = tatami::make_DelayedUnaryIsometricOp(input->ptr, tatami::DelayedCosHelper<>());
    } else if (op == "sin") {
        output->ptr = tatami::make_DelayedUnaryIsometricOp(input->ptr, tatami::DelayedSinHelper<>());
    } else if (op == "tan") {
        output->ptr = tatami::make_DelayedUnaryIsometricOp(input->ptr, tatami::DelayedTanHelper<>());
    } else if (op == "cospi") {
        output->ptr = tatami::make_DelayedUnaryIsometricOp(
            tatami::make_DelayedUnaryIsometricOp(input->ptr, tatami::make_DelayedMultiplyScalarHelper(M_PI)),
            tatami::DelayedCosHelper<>()
        );
    } else if (op == "sinpi") {
        output->ptr = tatami::make_DelayedUnaryIsometricOp(
            tatami::make_DelayedUnaryIsometricOp(input->ptr, tatami::make_DelayedMultiplyScalarHelper(M_PI)),
            tatami::DelayedSinHelper<>()
        );
    } else if (op == "tanpi") {
        output->ptr = tatami::make_DelayedUnaryIsometricOp(
            tatami::make_DelayedUnaryIsometricOp(input->ptr, tatami::make_DelayedMultiplyScalarHelper(M_PI)),
            tatami::DelayedTanHelper<>()
        );
    } else if (op == "acos") {
        output->ptr = tatami::make_DelayedUnaryIsometricOp(input->ptr, tatami::DelayedAcosHelper<>());
    } else if (op == "asin") {
        output->ptr = tatami::make_DelayedUnaryIsometricOp(input->ptr, tatami::DelayedAsinHelper<>());
    } else if (op == "atan") {
        output->ptr = tatami::make_DelayedUnaryIsometricOp(input->ptr, tatami::DelayedAtanHelper<>());
    } else if (op == "cosh") {
        output->ptr = tatami::make_DelayedUnaryIsometricOp(input->ptr, tatami::DelayedCoshHelper<>());
    } else if (op == "sinh") {
        output->ptr = tatami::make_DelayedUnaryIsometricOp(input->ptr, tatami::DelayedSinhHelper<>());
    } else if (op == "tanh") {
        output->ptr = tatami::make_DelayedUnaryIsometricOp(input->ptr, tatami::DelayedTanhHelper<>());
    } else if (op == "acosh") {
        output->ptr = tatami::make_DelayedUnaryIsometricOp(input->ptr, tatami::DelayedAcoshHelper<>());
    } else if (op == "asinh") {
        output->ptr = tatami::make_DelayedUnaryIsometricOp(input->ptr, tatami::DelayedAsinhHelper<>());
    } else if (op == "atanh") {
        output->ptr = tatami::make_DelayedUnaryIsometricOp(input->ptr, tatami::DelayedAtanhHelper<>());

    } else if (op == "lgamma") {
        output->ptr = tatami::make_DelayedUnaryIsometricOp(input->ptr, tatami::DelayedLgammaHelper<>());
    } else if (op == "gamma") {
        output->ptr = tatami::make_DelayedUnaryIsometricOp(input->ptr, tatami::DelayedGammaHelper<>());
    } else {
        return(R_NilValue);
    }

    output->original = input->original; // copying the reference to propagate GC protection.
    return output;
}

//[[Rcpp::export(rng=false)]]
SEXP apply_delayed_round(SEXP raw_input) {
    Rtatami::BoundNumericPointer input(raw_input);
    auto output = Rtatami::new_BoundNumericMatrix();
    output->ptr = tatami::make_DelayedUnaryIsometricOp(input->ptr, tatami::DelayedRoundHelper<>());
    output->original = input->original; // copying the reference to propagate GC protection.
    return output;
}
