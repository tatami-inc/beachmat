#include "Rtatami.h"
#include "Rcpp.h"
#include "tatami/tatami.hpp"
#include "Rmath.h"

//[[Rcpp::export(rng=false)]]
SEXP apply_delayed_log(SEXP raw_input, double base) {
    Rtatami::BoundNumericPointer input(raw_input);
    auto output = Rtatami::new_BoundNumericMatrix();
    output->ptr = tatami::make_DelayedUnaryIsometricOperation(input->ptr, tatami::DelayedUnaryIsometricLog(base));
    output->original = input->original; // copying the reference to propagate GC protection.
    return output;
}

//[[Rcpp::export(rng=false)]]
SEXP apply_delayed_unary_math(SEXP raw_input, const std::string& op) {
    Rtatami::BoundNumericPointer input(raw_input);
    auto output = Rtatami::new_BoundNumericMatrix();

    if (op == "abs") {
        output->ptr = tatami::make_DelayedUnaryIsometricOperation(input->ptr, tatami::DelayedUnaryIsometricAbs<>());
    } else if (op == "sign") {
        output->ptr = tatami::make_DelayedUnaryIsometricOperation(input->ptr, tatami::DelayedUnaryIsometricSign<>());
    } else if (op == "sqrt") {
        output->ptr = tatami::make_DelayedUnaryIsometricOperation(input->ptr, tatami::DelayedUnaryIsometricSqrt<>());
    } else if (op == "floor") {
        output->ptr = tatami::make_DelayedUnaryIsometricOperation(input->ptr, tatami::DelayedUnaryIsometricFloor<>());
    } else if (op == "ceiling") {
        output->ptr = tatami::make_DelayedUnaryIsometricOperation(input->ptr, tatami::DelayedUnaryIsometricCeiling<>());
    } else if (op == "trunc") {
        output->ptr = tatami::make_DelayedUnaryIsometricOperation(input->ptr, tatami::DelayedUnaryIsometricTrunc<>());

    } else if (op == "exp") {
        output->ptr = tatami::make_DelayedUnaryIsometricOperation(input->ptr, tatami::DelayedUnaryIsometricExp<>());
    } else if (op == "expm1") {
        output->ptr = tatami::make_DelayedUnaryIsometricOperation(input->ptr, tatami::DelayedUnaryIsometricExpm1<>());
    } else if (op == "log1p") {
        output->ptr = tatami::make_DelayedUnaryIsometricOperation(input->ptr, tatami::DelayedUnaryIsometricLog1p<>());
    } else if (op == "cos") {
        output->ptr = tatami::make_DelayedUnaryIsometricOperation(input->ptr, tatami::DelayedUnaryIsometricCos<>());
    } else if (op == "sin") {
        output->ptr = tatami::make_DelayedUnaryIsometricOperation(input->ptr, tatami::DelayedUnaryIsometricSin<>());
    } else if (op == "tan") {
        output->ptr = tatami::make_DelayedUnaryIsometricOperation(input->ptr, tatami::DelayedUnaryIsometricTan<>());
    } else if (op == "cospi") {
        output->ptr = tatami::make_DelayedUnaryIsometricOperation(
            tatami::make_DelayedUnaryIsometricOperation(input->ptr, tatami::make_DelayedUnaryIsometricMultiplyScalar(M_PI)),
            tatami::DelayedUnaryIsometricCos<>()
        );
    } else if (op == "sinpi") {
        output->ptr = tatami::make_DelayedUnaryIsometricOperation(
            tatami::make_DelayedUnaryIsometricOperation(input->ptr, tatami::make_DelayedUnaryIsometricMultiplyScalar(M_PI)),
            tatami::DelayedUnaryIsometricSin<>()
        );
    } else if (op == "tanpi") {
        output->ptr = tatami::make_DelayedUnaryIsometricOperation(
            tatami::make_DelayedUnaryIsometricOperation(input->ptr, tatami::make_DelayedUnaryIsometricMultiplyScalar(M_PI)),
            tatami::DelayedUnaryIsometricTan<>()
        );
    } else if (op == "acos") {
        output->ptr = tatami::make_DelayedUnaryIsometricOperation(input->ptr, tatami::DelayedUnaryIsometricAcos<>());
    } else if (op == "asin") {
        output->ptr = tatami::make_DelayedUnaryIsometricOperation(input->ptr, tatami::DelayedUnaryIsometricAsin<>());
    } else if (op == "atan") {
        output->ptr = tatami::make_DelayedUnaryIsometricOperation(input->ptr, tatami::DelayedUnaryIsometricAtan<>());
    } else if (op == "cosh") {
        output->ptr = tatami::make_DelayedUnaryIsometricOperation(input->ptr, tatami::DelayedUnaryIsometricCosh<>());
    } else if (op == "sinh") {
        output->ptr = tatami::make_DelayedUnaryIsometricOperation(input->ptr, tatami::DelayedUnaryIsometricSinh<>());
    } else if (op == "tanh") {
        output->ptr = tatami::make_DelayedUnaryIsometricOperation(input->ptr, tatami::DelayedUnaryIsometricTanh<>());
    } else if (op == "acosh") {
        output->ptr = tatami::make_DelayedUnaryIsometricOperation(input->ptr, tatami::DelayedUnaryIsometricAcosh<>());
    } else if (op == "asinh") {
        output->ptr = tatami::make_DelayedUnaryIsometricOperation(input->ptr, tatami::DelayedUnaryIsometricAsinh<>());
    } else if (op == "atanh") {
        output->ptr = tatami::make_DelayedUnaryIsometricOperation(input->ptr, tatami::DelayedUnaryIsometricAtanh<>());

    } else if (op == "lgamma") {
        output->ptr = tatami::make_DelayedUnaryIsometricOperation(input->ptr, tatami::DelayedUnaryIsometricLgamma<>());
    } else if (op == "gamma") {
        output->ptr = tatami::make_DelayedUnaryIsometricOperation(input->ptr, tatami::DelayedUnaryIsometricGamma<>());
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
    output->ptr = tatami::make_DelayedUnaryIsometricOperation(input->ptr, tatami::DelayedUnaryIsometricRound<>());
    output->original = input->original; // copying the reference to propagate GC protection.
    return output;
}
