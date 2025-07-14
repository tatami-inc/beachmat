#include "Rtatami.h"
#include "Rcpp.h"
#include "tatami/tatami.hpp"
#include "Rmath.h"

#include <memory>
#include <string>
#include <stdexcept>

//[[Rcpp::export(rng=false)]]
SEXP apply_delayed_associative_arithmetic(SEXP raw_input, Rcpp::NumericVector val, bool row, std::string op) {
    Rtatami::BoundNumericPointer input(raw_input);
    Rcpp::List protectorate(2);
    protectorate[0] = input->original;

    std::shared_ptr<tatami::DelayedUnaryIsometricOperationHelper<double, double, int> > opptr;
    if (val.size() == 1) {
        protectorate[1] = R_NilValue;
        if (op == "+") {
            opptr.reset(new tatami::DelayedUnaryIsometricAddScalarHelper<double, double, int, double>(val[0]));
        } else if (op == "*") {
            opptr.reset(new tatami::DelayedUnaryIsometricMultiplyScalarHelper<double, double, int, double>(val[0]));
        } else {
            throw std::runtime_error("unknown associative arithmetic operation '" + op + "'");
        }

    } else {
        protectorate[1] = val;
        tatami::ArrayView<double> view(static_cast<const double*>(val.begin()), val.size());
        if (op == "+") {
            opptr.reset(new tatami::DelayedUnaryIsometricAddVectorHelper<double, double, int, decltype(view)>(std::move(view), row));
        } else if (op == "*") {
            opptr.reset(new tatami::DelayedUnaryIsometricMultiplyVectorHelper<double, double, int, decltype(view)>(std::move(view), row));
        } else {
            throw std::runtime_error("unknown associative arithmetic operation '" + op + "'");
        }
    }

    auto output = Rtatami::new_BoundNumericMatrix();
    output->ptr.reset(new tatami::DelayedUnaryIsometricOperation<double, double, int>(input->ptr, std::move(opptr)));
    output->original = protectorate; // propagate protection for all child objects by copying references.
    return output;
}

template<bool right_>
std::shared_ptr<tatami::DelayedUnaryIsometricOperationHelper<double, double, int> > set_delayed_nonassociative_arithmetic_scalar(
    double val,
    const std::string& op)
{
    std::shared_ptr<tatami::DelayedUnaryIsometricOperationHelper<double, double, int> > opptr;
    if (op == "-") {
        opptr.reset(new tatami::DelayedUnaryIsometricSubtractScalarHelper<right_, double, double, int, double>(val));
    } else if (op == "/") {
        opptr.reset(new tatami::DelayedUnaryIsometricDivideScalarHelper<right_, double, double, int, double>(val));
    } else if (op == "%/%") {
        opptr.reset(new tatami::DelayedUnaryIsometricIntegerDivideScalarHelper<right_, double, double, int, double>(val));
    } else if (op == "^") {
        opptr.reset(new tatami::DelayedUnaryIsometricPowerScalarHelper<right_, double, double, int, double>(val));
    } else if (op == "%%") {
        opptr.reset(new tatami::DelayedUnaryIsometricModuloScalarHelper<right_, double, double, int, double>(val));
    } else {
        throw std::runtime_error("unknown non-associative arithmetic operation '" + op + "'");
    }
    return opptr;
}

template<bool right_>
std::shared_ptr<tatami::DelayedUnaryIsometricOperationHelper<double, double, int> > set_delayed_nonassociative_arithmetic_vector(
    const Rcpp::NumericVector& val,
    const std::string& op,
    bool row)
{
    std::shared_ptr<tatami::DelayedUnaryIsometricOperationHelper<double, double, int> > opptr;
    tatami::ArrayView<double> view(static_cast<const double*>(val.begin()), val.size());
    if (op == "-") {
        opptr.reset(new tatami::DelayedUnaryIsometricSubtractVectorHelper<right_, double, double, int, decltype(view)>(std::move(view), row));
    } else if (op == "/") {
        opptr.reset(new tatami::DelayedUnaryIsometricDivideVectorHelper<right_, double, double, int, decltype(view)>(std::move(view), row));
    } else if (op == "%/%") {
        opptr.reset(new tatami::DelayedUnaryIsometricIntegerDivideVectorHelper<right_, double, double, int, decltype(view)>(std::move(view), row));
    } else if (op == "^") {
        opptr.reset(new tatami::DelayedUnaryIsometricPowerVectorHelper<right_, double, double, int, decltype(view)>(std::move(view), row));
    } else if (op == "%%") {
        opptr.reset(new tatami::DelayedUnaryIsometricModuloVectorHelper<right_, double, double, int, decltype(view)>(std::move(view), row));
    } else {
        throw std::runtime_error("unknown non-associative arithmetic operation '" + op + "'");
    }
    return opptr;
}

//[[Rcpp::export(rng=false)]]
SEXP apply_delayed_nonassociative_arithmetic(SEXP raw_input, Rcpp::NumericVector val, bool right, bool row, std::string op) {
    Rtatami::BoundNumericPointer input(raw_input);
    Rcpp::List protectorate(2);
    protectorate[0] = input->original;

    std::shared_ptr<tatami::DelayedUnaryIsometricOperationHelper<double, double, int> > opptr;
    if (val.size() == 1) {
        protectorate[1] = R_NilValue;
        if (right) {
            opptr = set_delayed_nonassociative_arithmetic_scalar<true>(val[0], op);
        } else {
            opptr = set_delayed_nonassociative_arithmetic_scalar<false>(val[0], op);
        }
    } else {
        protectorate[1] = val;        
        if (right) {
            opptr = set_delayed_nonassociative_arithmetic_vector<true>(val, op, row);
        } else {
            opptr = set_delayed_nonassociative_arithmetic_vector<false>(val, op, row);
        }
    }

    auto output = Rtatami::new_BoundNumericMatrix();
    output->ptr.reset(new tatami::DelayedUnaryIsometricOperation<double, double, int>(input->ptr, std::move(opptr)));
    output->original = protectorate; // propagate protection for all child objects by copying references.
    return output;
}

//[[Rcpp::export(rng=false)]]
SEXP apply_delayed_comparison(SEXP raw_input, Rcpp::NumericVector val, bool row, std::string op) {
    Rtatami::BoundNumericPointer input(raw_input);
    const auto& shared = input->ptr;
    Rcpp::List protectorate(2);
    protectorate[0] = input->original;

    std::shared_ptr<tatami::DelayedUnaryIsometricOperationHelper<double, double, int> > opptr;
    if (val.size() == 1) {
        protectorate[1] = R_NilValue;
        if (op == "==") {
            opptr.reset(new tatami::DelayedUnaryIsometricEqualScalarHelper<double, double, int, double>(val[0]));
        } else if (op == ">") {
            opptr.reset(new tatami::DelayedUnaryIsometricGreaterThanScalarHelper<double, double, int, double>(val[0]));
        } else if (op == "<") {
            opptr.reset(new tatami::DelayedUnaryIsometricLessThanScalarHelper<double, double, int, double>(val[0]));
        } else if (op == ">=") {
            opptr.reset(new tatami::DelayedUnaryIsometricGreaterThanOrEqualScalarHelper<double, double, int, double>(val[0]));
        } else if (op == "<=") {
            opptr.reset(new tatami::DelayedUnaryIsometricLessThanOrEqualScalarHelper<double, double, int, double>(val[0]));
        } else if (op == "!=") {
            opptr.reset(new tatami::DelayedUnaryIsometricNotEqualScalarHelper<double, double, int, double>(val[0]));
        } else {
            throw std::runtime_error("unknown delayed comparison operation '" + op + "'");
        }
    } else {
        protectorate[1] = val;
        tatami::ArrayView<double> view(static_cast<const double*>(val.begin()), val.size());
        if (op == "==") {
            opptr.reset(new tatami::DelayedUnaryIsometricEqualVectorHelper<double, double, int, decltype(view)>(std::move(view), row));
        } else if (op == ">") {
            opptr.reset(new tatami::DelayedUnaryIsometricGreaterThanVectorHelper<double, double, int, decltype(view)>(std::move(view), row));
        } else if (op == "<") {
            opptr.reset(new tatami::DelayedUnaryIsometricLessThanVectorHelper<double, double, int, decltype(view)>(std::move(view), row));
        } else if (op == ">=") {
            opptr.reset(new tatami::DelayedUnaryIsometricGreaterThanOrEqualVectorHelper<double, double, int, decltype(view)>(std::move(view), row));
        } else if (op == "<=") {
            opptr.reset(new tatami::DelayedUnaryIsometricLessThanOrEqualVectorHelper<double, double, int, decltype(view)>(std::move(view), row));
        } else if (op == "!=") {
            opptr.reset(new tatami::DelayedUnaryIsometricNotEqualVectorHelper<double, double, int, decltype(view)>(std::move(view), row));
        } else {
            throw std::runtime_error("unknown delayed comparison operation '" + op + "'");
        }
    }

    auto output = Rtatami::new_BoundNumericMatrix();
    output->ptr.reset(new tatami::DelayedUnaryIsometricOperation<double, double, int>(input->ptr, std::move(opptr)));
    output->original = protectorate; // propagate protection for all child objects by copying references.
    return output;
}

//[[Rcpp::export(rng=false)]]
SEXP apply_delayed_boolean(SEXP raw_input, Rcpp::LogicalVector val, bool row, std::string op) {
    Rtatami::BoundNumericPointer input(raw_input);
    const auto& shared = input->ptr;
    Rcpp::List protectorate(2);
    protectorate[0] = input->original;

    std::shared_ptr<tatami::DelayedUnaryIsometricOperationHelper<double, double, int> > opptr;
    if (val.size() == 1) {
        protectorate[1] = R_NilValue;
        if (op == "&") {
            opptr.reset(new tatami::DelayedUnaryIsometricBooleanAndScalarHelper<double, double, int>(val[0]));
        } else if (op == "|") {
            opptr.reset(new tatami::DelayedUnaryIsometricBooleanOrScalarHelper<double, double, int>(val[0]));
        } else {
            throw std::runtime_error("unknown delayed boolean operation '" + op + "'");
        }
    } else {
        protectorate[1] = val;
        tatami::ArrayView<int> view(static_cast<const int*>(val.begin()), val.size());
        if (op == "&") {
            opptr.reset(new tatami::DelayedUnaryIsometricBooleanAndVectorHelper<double, double, int, decltype(view)>(std::move(view), row));
        } else if (op == "|") {
            opptr.reset(new tatami::DelayedUnaryIsometricBooleanOrVectorHelper<double, double, int, decltype(view)>(std::move(view), row));
        } else {
            throw std::runtime_error("unknown delayed boolean operation '" + op + "'");
        }
    }

    auto output = Rtatami::new_BoundNumericMatrix();
    output->ptr.reset(new tatami::DelayedUnaryIsometricOperation<double, double, int>(input->ptr, std::move(opptr)));
    output->original = protectorate; // propagate protection for all child objects by copying references.
    return output;
}

//[[Rcpp::export(rng=false)]]
SEXP apply_delayed_boolean_not(SEXP raw_input) {
    Rtatami::BoundNumericPointer input(raw_input);
    auto output = Rtatami::new_BoundNumericMatrix();
    output->ptr.reset(new tatami::DelayedUnaryIsometricOperation<double, double, int>(input->ptr, std::make_shared<tatami::DelayedUnaryIsometricBooleanNotHelper<double, double, int> >()));
    output->original = input->original; // copying the reference to propagate GC protection.
    return output;
}
