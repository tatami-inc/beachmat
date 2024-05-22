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
        output->ptr = tatami::make_DelayedBinaryIsometricOperation(left_shared, right_shared, tatami::make_DelayedBinaryIsometricAdd());

    } else if (op == "-") {
        output->ptr = tatami::make_DelayedBinaryIsometricOperation(left_shared, right_shared, tatami::make_DelayedBinaryIsometricSubtract());

    } else if (op == "*") {
        output->ptr = tatami::make_DelayedBinaryIsometricOperation(left_shared, right_shared, tatami::make_DelayedBinaryIsometricMultiply());

    } else if (op == "/") {
        output->ptr = tatami::make_DelayedBinaryIsometricOperation(left_shared, right_shared, tatami::make_DelayedBinaryIsometricDivide());

    } else if (op == "%/%") {
        output->ptr = tatami::make_DelayedBinaryIsometricOperation(left_shared, right_shared, tatami::make_DelayedBinaryIsometricIntegerDivide());

    } else if (op == "^") {
        output->ptr = tatami::make_DelayedBinaryIsometricOperation(left_shared, right_shared, tatami::make_DelayedBinaryIsometricPower());

    } else if (op == "%%") {
        output->ptr = tatami::make_DelayedBinaryIsometricOperation(left_shared, right_shared, tatami::make_DelayedBinaryIsometricModulo());

    } else if (op == "==") {
        output->ptr = tatami::make_DelayedBinaryIsometricOperation(left_shared, right_shared, tatami::make_DelayedBinaryIsometricEqual());

    } else if (op == "!=") {
        output->ptr = tatami::make_DelayedBinaryIsometricOperation(left_shared, right_shared, tatami::make_DelayedBinaryIsometricNotEqual());

    } else if (op == ">") {
        output->ptr = tatami::make_DelayedBinaryIsometricOperation(left_shared, right_shared, tatami::make_DelayedBinaryIsometricGreaterThan());

    } else if (op == "<") {
        output->ptr = tatami::make_DelayedBinaryIsometricOperation(left_shared, right_shared, tatami::make_DelayedBinaryIsometricLessThan());

    } else if (op == ">=") {
        output->ptr = tatami::make_DelayedBinaryIsometricOperation(left_shared, right_shared, tatami::make_DelayedBinaryIsometricGreaterThanOrEqual());

    } else if (op == "<=") {
        output->ptr = tatami::make_DelayedBinaryIsometricOperation(left_shared, right_shared, tatami::make_DelayedBinaryIsometricLessThanOrEqual());

    } else if (op == "&") {
        output->ptr = tatami::make_DelayedBinaryIsometricOperation(left_shared, right_shared, tatami::make_DelayedBinaryIsometricBooleanAnd());

    } else if (op == "|") {
        output->ptr = tatami::make_DelayedBinaryIsometricOperation(left_shared, right_shared, tatami::make_DelayedBinaryIsometricBooleanOr());

    } else {
        throw std::runtime_error("unknown delayed binary operation '" + op + "'");
    }

    output->original = protectorate; // propagate protection for all child objects by copying references.
    return output;
}
