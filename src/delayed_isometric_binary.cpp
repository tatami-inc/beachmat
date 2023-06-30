#include "Rtatami.h"
#include "Rcpp.h"
#include "tatami/tatami.hpp"

//[[Rcpp::export(rng=false)]]
SEXP apply_delayed_binary_operation(SEXP left_input, SEXP right_input, std::string op) {
    Rtatami::BoundNumericPointer left(left_input);
    Rtatami::BoundNumericPointer right(right_input);
    const auto& left_shared = left->ptr;
    const auto& right_shared = right->ptr;
    Rcpp::List protectorate(2);
    protectorate[0] = left->original;
    protectorate[1] = right->original;

    auto output = Rtatami::new_BoundNumericMatrix();
    if (op == "+") {
        output->ptr = tatami::make_DelayedBinaryIsometricOp(left_shared, right_shared, tatami::make_DelayedBinaryAddHelper());

    } else if (op == "-") {
        output->ptr = tatami::make_DelayedBinaryIsometricOp(left_shared, right_shared, tatami::make_DelayedBinarySubtractHelper());

    } else if (op == "*") {
        output->ptr = tatami::make_DelayedBinaryIsometricOp(left_shared, right_shared, tatami::make_DelayedBinaryMultiplyHelper());

    } else if (op == "/") {
        output->ptr = tatami::make_DelayedBinaryIsometricOp(left_shared, right_shared, tatami::make_DelayedBinaryDivideHelper());

    } else if (op == "%/%") {
        output->ptr = tatami::make_DelayedBinaryIsometricOp(left_shared, right_shared, tatami::make_DelayedBinaryIntegerDivideHelper());

    } else if (op == "^") {
        output->ptr = tatami::make_DelayedBinaryIsometricOp(left_shared, right_shared, tatami::make_DelayedBinaryPowerHelper());

    } else if (op == "%%") {
        output->ptr = tatami::make_DelayedBinaryIsometricOp(left_shared, right_shared, tatami::make_DelayedBinaryModuloHelper());

    } else if (op == "==") {
        output->ptr = tatami::make_DelayedBinaryIsometricOp(left_shared, right_shared, tatami::make_DelayedBinaryEqualHelper());

    } else if (op == "!=") {
        output->ptr = tatami::make_DelayedBinaryIsometricOp(left_shared, right_shared, tatami::make_DelayedBinaryNotEqualHelper());

    } else if (op == ">") {
        output->ptr = tatami::make_DelayedBinaryIsometricOp(left_shared, right_shared, tatami::make_DelayedBinaryGreaterThanHelper());

    } else if (op == "<") {
        output->ptr = tatami::make_DelayedBinaryIsometricOp(left_shared, right_shared, tatami::make_DelayedBinaryLessThanHelper());

    } else if (op == ">=") {
        output->ptr = tatami::make_DelayedBinaryIsometricOp(left_shared, right_shared, tatami::make_DelayedBinaryGreaterThanOrEqualHelper());

    } else if (op == "<=") {
        output->ptr = tatami::make_DelayedBinaryIsometricOp(left_shared, right_shared, tatami::make_DelayedBinaryLessThanOrEqualHelper());

    } else if (op == "&") {
        output->ptr = tatami::make_DelayedBinaryIsometricOp(left_shared, right_shared, tatami::make_DelayedBinaryBooleanAndHelper());

    } else if (op == "|") {
        output->ptr = tatami::make_DelayedBinaryIsometricOp(left_shared, right_shared, tatami::make_DelayedBinaryBooleanOrHelper());

    } else {
        throw std::runtime_error("unknown delayed binary operation '" + op + "'");
    }

    output->original = protectorate; // propagate protection for all child objects by copying references.
    return output;
}
