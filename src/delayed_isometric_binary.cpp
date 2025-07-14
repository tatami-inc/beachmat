#include "Rtatami.h"
#include "Rcpp.h"
#include "tatami/tatami.hpp"

#include <memory>
#include <string>
#include <stdexcept>

//[[Rcpp::export(rng=false)]]
SEXP apply_delayed_binary_operation(SEXP left_input, SEXP right_input, std::string op) {
    Rtatami::BoundNumericPointer left(left_input);
    Rtatami::BoundNumericPointer right(right_input);
    const auto& left_shared = left->ptr;
    const auto& right_shared = right->ptr;
    Rcpp::List protectorate(2);
    protectorate[0] = left->original;
    protectorate[1] = right->original;

    std::shared_ptr<tatami::DelayedBinaryIsometricOperationHelper<double, double, int> > opptr;
    if (op == "+") {
        opptr.reset(new tatami::DelayedBinaryIsometricAddHelper<double, double, int>);
    } else if (op == "-") {
        opptr.reset(new tatami::DelayedBinaryIsometricSubtractHelper<double, double, int>);
    } else if (op == "*") {
        opptr.reset(new tatami::DelayedBinaryIsometricMultiplyHelper<double, double, int>);
    } else if (op == "/") {
        opptr.reset(new tatami::DelayedBinaryIsometricDivideHelper<double, double, int>);
    } else if (op == "%/%") {
        opptr.reset(new tatami::DelayedBinaryIsometricIntegerDivideHelper<double, double, int>);
    } else if (op == "^") {
        opptr.reset(new tatami::DelayedBinaryIsometricPowerHelper<double, double, int>);
    } else if (op == "%%") {
        opptr.reset(new tatami::DelayedBinaryIsometricModuloHelper<double, double, int>);
    } else if (op == "==") {
        opptr.reset(new tatami::DelayedBinaryIsometricEqualHelper<double, double, int>);
    } else if (op == "!=") {
        opptr.reset(new tatami::DelayedBinaryIsometricNotEqualHelper<double, double, int>);
    } else if (op == ">") {
        opptr.reset(new tatami::DelayedBinaryIsometricGreaterThanHelper<double, double, int>);
    } else if (op == "<") {
        opptr.reset(new tatami::DelayedBinaryIsometricLessThanHelper<double, double, int>);
    } else if (op == ">=") {
        opptr.reset(new tatami::DelayedBinaryIsometricGreaterThanOrEqualHelper<double, double, int>);
    } else if (op == "<=") {
        opptr.reset(new tatami::DelayedBinaryIsometricLessThanOrEqualHelper<double, double, int>);
    } else if (op == "&") {
        opptr.reset(new tatami::DelayedBinaryIsometricBooleanAndHelper<double, double, int>);
    } else if (op == "|") {
        opptr.reset(new tatami::DelayedBinaryIsometricBooleanOrHelper<double, double, int>);
    } else {
        throw std::runtime_error("unknown delayed binary operation '" + op + "'");
    }

    auto output = Rtatami::new_BoundNumericMatrix();
    output->ptr.reset(new tatami::DelayedBinaryIsometricOperation<double, double, int>(left_shared, right_shared, std::move(opptr)));
    output->original = protectorate; // propagate protection for all child objects by copying references.
    return output;
}
