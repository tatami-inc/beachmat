#include "Rtatami.h"
#include "Rcpp.h"
#include "tatami/tatami.hpp"
#include "Rmath.h"

void set_delayed_associative_arithmetic_vector(const std::shared_ptr<tatami::NumericMatrix>& shared, const Rcpp::NumericVector& val, const std::string& op, std::shared_ptr<tatami::NumericMatrix>& outptr, bool row) {
    tatami::ArrayView<double> view(static_cast<const double*>(val.begin()), val.size());
    if (op == "+") {
        outptr = tatami::make_DelayedUnaryIsometricOperation(shared, tatami::make_DelayedUnaryIsometricAddVector(std::move(view), row));
    } else if (op == "*") {
        outptr = tatami::make_DelayedUnaryIsometricOperation(shared, tatami::make_DelayedUnaryIsometricMultiplyVector(std::move(view), row));
    } else {
        throw std::runtime_error("unknown associative arithmetic operation '" + op + "'");
    }
}

//[[Rcpp::export(rng=false)]]
SEXP apply_delayed_associative_arithmetic(SEXP raw_input, Rcpp::NumericVector val, bool row, std::string op) {
    Rtatami::BoundNumericPointer input(raw_input);
    const auto& shared = input->ptr;
    Rcpp::List protectorate(2);
    protectorate[0] = input->original;

    auto output = Rtatami::new_BoundNumericMatrix();
    if (val.size() == 1) {
        if (op == "+") {
            output->ptr = tatami::make_DelayedUnaryIsometricOperation(shared, tatami::make_DelayedUnaryIsometricAddScalar(val[0]));
        } else if (op == "*") {
            output->ptr = tatami::make_DelayedUnaryIsometricOperation(shared, tatami::make_DelayedUnaryIsometricMultiplyScalar(val[0]));
        } else {
            throw std::runtime_error("unknown associative arithmetic operation '" + op + "'");
        }

    } else {
        protectorate[1] = val;
        set_delayed_associative_arithmetic_vector(shared, val, op, output->ptr, row);
    }

    output->original = protectorate; // propagate protection for all child objects by copying references.
    return output;
}

template<bool right_>
void set_delayed_nonassociative_arithmetic_scalar(const std::shared_ptr<tatami::NumericMatrix>& shared, double val, const std::string& op, std::shared_ptr<tatami::NumericMatrix>& outptr) {
    if (op == "-") {
        outptr = tatami::make_DelayedUnaryIsometricOperation(shared, tatami::make_DelayedUnaryIsometricSubtractScalar<right_>(val));
    } else if (op == "/") {
        outptr = tatami::make_DelayedUnaryIsometricOperation(shared, tatami::make_DelayedUnaryIsometricDivideScalar<right_>(val));
    } else if (op == "%/%") {
        outptr = tatami::make_DelayedUnaryIsometricOperation(shared, tatami::make_DelayedUnaryIsometricIntegerDivideScalar<right_>(val));
    } else if (op == "^") {
        outptr = tatami::make_DelayedUnaryIsometricOperation(shared, tatami::make_DelayedUnaryIsometricPowerScalar<right_>(val));
    } else if (op == "%%") {
        outptr = tatami::make_DelayedUnaryIsometricOperation(shared, tatami::make_DelayedUnaryIsometricModuloScalar<right_>(val));
    } else {
        throw std::runtime_error("unknown non-associative arithmetic operation '" + op + "'");
    }
}

template<bool right_>
void set_delayed_nonassociative_arithmetic_vector(const std::shared_ptr<tatami::NumericMatrix>& shared, const Rcpp::NumericVector& val, const std::string& op, std::shared_ptr<tatami::NumericMatrix>& outptr, bool row) {
    tatami::ArrayView<double> view(static_cast<const double*>(val.begin()), val.size());
    if (op == "-") {
        outptr = tatami::make_DelayedUnaryIsometricOperation(shared, tatami::make_DelayedUnaryIsometricSubtractVector<right_>(std::move(view), row));
    } else if (op == "/") {
        outptr = tatami::make_DelayedUnaryIsometricOperation(shared, tatami::make_DelayedUnaryIsometricDivideVector<right_>(std::move(view), row));
    } else if (op == "%/%") {
        outptr = tatami::make_DelayedUnaryIsometricOperation(shared, tatami::make_DelayedUnaryIsometricIntegerDivideVector<right_>(std::move(view), row));
    } else if (op == "^") {
        outptr = tatami::make_DelayedUnaryIsometricOperation(shared, tatami::make_DelayedUnaryIsometricPowerVector<right_>(std::move(view), row));
    } else if (op == "%%") {
        outptr = tatami::make_DelayedUnaryIsometricOperation(shared, tatami::make_DelayedUnaryIsometricModuloVector<right_>(std::move(view), row));
    } else {
        throw std::runtime_error("unknown non-associative arithmetic operation '" + op + "'");
    }
}

//[[Rcpp::export(rng=false)]]
SEXP apply_delayed_nonassociative_arithmetic(SEXP raw_input, Rcpp::NumericVector val, bool right, bool row, std::string op) {
    Rtatami::BoundNumericPointer input(raw_input);
    const auto& shared = input->ptr;
    Rcpp::List protectorate(2);
    protectorate[0] = input->original;

    auto output = Rtatami::new_BoundNumericMatrix();
    if (val.size() == 1) {
        if (right) {
            set_delayed_nonassociative_arithmetic_scalar<true>(shared, val[0], op, output->ptr);
        } else {
            set_delayed_nonassociative_arithmetic_scalar<false>(shared, val[0], op, output->ptr);
        }

    } else {
        protectorate[1] = val;        
        if (right) {
            set_delayed_nonassociative_arithmetic_vector<true>(shared, val, op, output->ptr, row);
        } else {
            set_delayed_nonassociative_arithmetic_vector<false>(shared, val, op, output->ptr, row);
        }
    }

    output->original = protectorate; // propagate protection for all child objects by copying references.
    return output;
}

void set_delayed_comparison_vector(const std::shared_ptr<tatami::NumericMatrix>& shared, const Rcpp::NumericVector& val, const std::string& op, std::shared_ptr<tatami::NumericMatrix>& outptr, bool row) {
    tatami::ArrayView<double> view(static_cast<const double*>(val.begin()), val.size());
    if (op == "==") {
        outptr = tatami::make_DelayedUnaryIsometricOperation(shared, tatami::make_DelayedUnaryIsometricEqualVector(std::move(view), row));
    } else if (op == ">" ){
        outptr = tatami::make_DelayedUnaryIsometricOperation(shared, tatami::make_DelayedUnaryIsometricGreaterThanVector(std::move(view), row));
    } else if (op == "<" ){
        outptr = tatami::make_DelayedUnaryIsometricOperation(shared, tatami::make_DelayedUnaryIsometricLessThanVector(std::move(view), row));
    } else if (op == ">=" ){
        outptr = tatami::make_DelayedUnaryIsometricOperation(shared, tatami::make_DelayedUnaryIsometricGreaterThanOrEqualVector(std::move(view), row));
    } else if (op == "<=" ){
        outptr = tatami::make_DelayedUnaryIsometricOperation(shared, tatami::make_DelayedUnaryIsometricLessThanOrEqualVector(std::move(view), row));
    } else if (op == "!=" ){
        outptr = tatami::make_DelayedUnaryIsometricOperation(shared, tatami::make_DelayedUnaryIsometricNotEqualVector(std::move(view), row));
    } else {
        throw std::runtime_error("unknown delayed comparison operation '" + op + "'");
    }
}

//[[Rcpp::export(rng=false)]]
SEXP apply_delayed_comparison(SEXP raw_input, Rcpp::NumericVector val, bool row, std::string op) {
    Rtatami::BoundNumericPointer input(raw_input);
    const auto& shared = input->ptr;
    Rcpp::List protectorate(2);
    protectorate[0] = input->original;

    auto output = Rtatami::new_BoundNumericMatrix();
    if (val.size() == 1) {
        if (op == "==") {
            output->ptr = tatami::make_DelayedUnaryIsometricOperation(shared, tatami::make_DelayedUnaryIsometricEqualScalar(val[0]));
        } else if (op == ">" ){
            output->ptr = tatami::make_DelayedUnaryIsometricOperation(shared, tatami::make_DelayedUnaryIsometricGreaterThanScalar(val[0]));
        } else if (op == "<" ){
            output->ptr = tatami::make_DelayedUnaryIsometricOperation(shared, tatami::make_DelayedUnaryIsometricLessThanScalar(val[0]));
        } else if (op == ">=" ){
            output->ptr = tatami::make_DelayedUnaryIsometricOperation(shared, tatami::make_DelayedUnaryIsometricGreaterThanOrEqualScalar(val[0]));
        } else if (op == "<=" ){
            output->ptr = tatami::make_DelayedUnaryIsometricOperation(shared, tatami::make_DelayedUnaryIsometricLessThanOrEqualScalar(val[0]));
        } else if (op == "!=" ){
            output->ptr = tatami::make_DelayedUnaryIsometricOperation(shared, tatami::make_DelayedUnaryIsometricNotEqualScalar(val[0]));
        } else {
            throw std::runtime_error("unknown delayed comparison operation '" + op + "'");
        }

    } else {
        protectorate[1] = val;
        set_delayed_comparison_vector(shared, val, op, output->ptr, row);
    }

    output->original = protectorate; // propagate protection for all child objects by copying references.
    return output;
}

void set_delayed_boolean_vector(const std::shared_ptr<tatami::NumericMatrix>& shared, const Rcpp::LogicalVector& val, const std::string& op, std::shared_ptr<tatami::NumericMatrix>& outptr, bool row) {
    tatami::ArrayView<int> view(static_cast<const int*>(val.begin()), val.size());
    if (op == "&") {
        outptr = tatami::make_DelayedUnaryIsometricOperation(shared, tatami::make_DelayedUnaryIsometricBooleanAndVector(std::move(view), row));
    } else if (op == "|" ){
        outptr = tatami::make_DelayedUnaryIsometricOperation(shared, tatami::make_DelayedUnaryIsometricBooleanOrVector(std::move(view), row));
    } else {
        throw std::runtime_error("unknown delayed boolean operation '" + op + "'");
    }
}

//[[Rcpp::export(rng=false)]]
SEXP apply_delayed_boolean(SEXP raw_input, Rcpp::LogicalVector val, bool row, std::string op) {
    Rtatami::BoundNumericPointer input(raw_input);
    const auto& shared = input->ptr;
    Rcpp::List protectorate(2);
    protectorate[0] = input->original;

    auto output = Rtatami::new_BoundNumericMatrix();
    if (val.size() == 1) {
        if (op == "&") {
            output->ptr = tatami::make_DelayedUnaryIsometricOperation(shared, tatami::make_DelayedUnaryIsometricBooleanAndScalar(val[0]));
        } else if (op == "|" ){
            output->ptr = tatami::make_DelayedUnaryIsometricOperation(shared, tatami::make_DelayedUnaryIsometricBooleanOrScalar(val[0]));
        } else {
            throw std::runtime_error("unknown delayed boolean operation '" + op + "'");
        }

    } else {
        protectorate[1] = val;
        set_delayed_boolean_vector(shared, val, op, output->ptr, row);
    }

    output->original = protectorate; // propagate protection for all child objects by copying references.
    return output;
}

//[[Rcpp::export(rng=false)]]
SEXP apply_delayed_boolean_not(SEXP raw_input) {
    Rtatami::BoundNumericPointer input(raw_input);
    auto output = Rtatami::new_BoundNumericMatrix();
    output->ptr = tatami::make_DelayedUnaryIsometricOperation(input->ptr, tatami::make_DelayedUnaryIsometricBooleanNot());
    output->original = input->original; // copying the reference to propagate GC protection.
    return output;
}
