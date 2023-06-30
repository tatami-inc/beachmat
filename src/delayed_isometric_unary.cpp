#include "Rtatami.h"
#include "Rcpp.h"
#include "tatami/tatami.hpp"
#include "Rmath.h"

template<int margin_>
void set_delayed_associative_arithmetic_vector(const std::shared_ptr<tatami::NumericMatrix>& shared, const Rcpp::NumericVector& val, const std::string& op, std::shared_ptr<tatami::NumericMatrix>& outptr) {
    tatami::ArrayView<double> view(static_cast<const double*>(val.begin()), val.size());
    if (op == "+") {
        outptr = tatami::make_DelayedUnaryIsometricOp(shared, tatami::make_DelayedAddVectorHelper<margin_>(std::move(view)));
    } else if (op == "*") {
        outptr = tatami::make_DelayedUnaryIsometricOp(shared, tatami::make_DelayedMultiplyVectorHelper<margin_>(std::move(view)));
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
            output->ptr = tatami::make_DelayedUnaryIsometricOp(shared, tatami::make_DelayedAddScalarHelper(val[0]));
        } else if (op == "*") {
            output->ptr = tatami::make_DelayedUnaryIsometricOp(shared, tatami::make_DelayedMultiplyScalarHelper(val[0]));
        } else {
            throw std::runtime_error("unknown associative arithmetic operation '" + op + "'");
        }

    } else {
        protectorate[1] = val;
        if (row) {
            set_delayed_associative_arithmetic_vector<0>(shared, val, op, output->ptr);
        } else {
            set_delayed_associative_arithmetic_vector<1>(shared, val, op, output->ptr);
        }
    }

    output->original = protectorate; // propagate protection for all child objects by copying references.
    return output;
}

template<bool right_>
void set_delayed_nonassociative_arithmetic_scalar(const std::shared_ptr<tatami::NumericMatrix>& shared, double val, const std::string& op, std::shared_ptr<tatami::NumericMatrix>& outptr) {
    if (op == "-") {
        outptr = tatami::make_DelayedUnaryIsometricOp(shared, tatami::make_DelayedSubtractScalarHelper<right_>(val));
    } else if (op == "/") {
        outptr = tatami::make_DelayedUnaryIsometricOp(shared, tatami::make_DelayedDivideScalarHelper<right_>(val));
    } else if (op == "%/%") {
        outptr = tatami::make_DelayedUnaryIsometricOp(shared, tatami::make_DelayedIntegerDivideScalarHelper<right_>(val));
    } else if (op == "^") {
        outptr = tatami::make_DelayedUnaryIsometricOp(shared, tatami::make_DelayedPowerScalarHelper<right_>(val));
    } else if (op == "%%") {
        outptr = tatami::make_DelayedUnaryIsometricOp(shared, tatami::make_DelayedModuloScalarHelper<right_>(val));
    } else {
        throw std::runtime_error("unknown non-associative arithmetic operation '" + op + "'");
    }
}

template<bool right_, int margin_>
void set_delayed_nonassociative_arithmetic_vector(const std::shared_ptr<tatami::NumericMatrix>& shared, const Rcpp::NumericVector& val, const std::string& op, std::shared_ptr<tatami::NumericMatrix>& outptr) {
    tatami::ArrayView<double> view(static_cast<const double*>(val.begin()), val.size());
    if (op == "-") {
        outptr = tatami::make_DelayedUnaryIsometricOp(shared, tatami::make_DelayedSubtractVectorHelper<right_, margin_>(std::move(view)));
    } else if (op == "/") {
        outptr = tatami::make_DelayedUnaryIsometricOp(shared, tatami::make_DelayedDivideVectorHelper<right_, margin_>(std::move(view)));
    } else if (op == "%/%") {
        outptr = tatami::make_DelayedUnaryIsometricOp(shared, tatami::make_DelayedIntegerDivideVectorHelper<right_, margin_>(std::move(view)));
    } else if (op == "^") {
        outptr = tatami::make_DelayedUnaryIsometricOp(shared, tatami::make_DelayedPowerVectorHelper<right_, margin_>(std::move(view)));
    } else if (op == "%%") {
        outptr = tatami::make_DelayedUnaryIsometricOp(shared, tatami::make_DelayedModuloVectorHelper<right_, margin_>(std::move(view)));
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
            if (row) {
                set_delayed_nonassociative_arithmetic_vector<true, 0>(shared, val, op, output->ptr);
            } else {
                set_delayed_nonassociative_arithmetic_vector<true, 1>(shared, val, op, output->ptr);
            }
        } else {
            if (row) {
                set_delayed_nonassociative_arithmetic_vector<false, 0>(shared, val, op, output->ptr);
            } else {
                set_delayed_nonassociative_arithmetic_vector<false, 1>(shared, val, op, output->ptr);
            }
        }
    }

    output->original = protectorate; // propagate protection for all child objects by copying references.
    return output;
}

template<int margin_>
void set_delayed_comparison_vector(const std::shared_ptr<tatami::NumericMatrix>& shared, const Rcpp::NumericVector& val, const std::string& op, std::shared_ptr<tatami::NumericMatrix>& outptr) {
    tatami::ArrayView<double> view(static_cast<const double*>(val.begin()), val.size());
    if (op == "==") {
        outptr = tatami::make_DelayedUnaryIsometricOp(shared, tatami::make_DelayedEqualVectorHelper<margin_>(std::move(view)));
    } else if (op == ">" ){
        outptr = tatami::make_DelayedUnaryIsometricOp(shared, tatami::make_DelayedGreaterThanVectorHelper<margin_>(std::move(view)));
    } else if (op == "<" ){
        outptr = tatami::make_DelayedUnaryIsometricOp(shared, tatami::make_DelayedLessThanVectorHelper<margin_>(std::move(view)));
    } else if (op == ">=" ){
        outptr = tatami::make_DelayedUnaryIsometricOp(shared, tatami::make_DelayedGreaterThanOrEqualVectorHelper<margin_>(std::move(view)));
    } else if (op == "<=" ){
        outptr = tatami::make_DelayedUnaryIsometricOp(shared, tatami::make_DelayedLessThanOrEqualVectorHelper<margin_>(std::move(view)));
    } else if (op == "!=" ){
        outptr = tatami::make_DelayedUnaryIsometricOp(shared, tatami::make_DelayedNotEqualVectorHelper<margin_>(std::move(view)));
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
            output->ptr = tatami::make_DelayedUnaryIsometricOp(shared, tatami::make_DelayedEqualScalarHelper(val[0]));
        } else if (op == ">" ){
            output->ptr = tatami::make_DelayedUnaryIsometricOp(shared, tatami::make_DelayedGreaterThanScalarHelper(val[0]));
        } else if (op == "<" ){
            output->ptr = tatami::make_DelayedUnaryIsometricOp(shared, tatami::make_DelayedLessThanScalarHelper(val[0]));
        } else if (op == ">=" ){
            output->ptr = tatami::make_DelayedUnaryIsometricOp(shared, tatami::make_DelayedGreaterThanOrEqualScalarHelper(val[0]));
        } else if (op == "<=" ){
            output->ptr = tatami::make_DelayedUnaryIsometricOp(shared, tatami::make_DelayedLessThanOrEqualScalarHelper(val[0]));
        } else if (op == "!=" ){
            output->ptr = tatami::make_DelayedUnaryIsometricOp(shared, tatami::make_DelayedNotEqualScalarHelper(val[0]));
        } else {
            throw std::runtime_error("unknown delayed comparison operation '" + op + "'");
        }

    } else {
        protectorate[1] = val;
        if (row) {
            set_delayed_comparison_vector<0>(shared, val, op, output->ptr);
        } else {
            set_delayed_comparison_vector<1>(shared, val, op, output->ptr);
        }
    }

    output->original = protectorate; // propagate protection for all child objects by copying references.
    return output;
}

template<int margin_>
void set_delayed_boolean_vector(const std::shared_ptr<tatami::NumericMatrix>& shared, const Rcpp::LogicalVector& val, const std::string& op, std::shared_ptr<tatami::NumericMatrix>& outptr) {
    tatami::ArrayView<int> view(static_cast<const int*>(val.begin()), val.size());
    if (op == "&") {
        outptr = tatami::make_DelayedUnaryIsometricOp(shared, tatami::make_DelayedBooleanAndVectorHelper<margin_>(std::move(view)));
    } else if (op == "|" ){
        outptr = tatami::make_DelayedUnaryIsometricOp(shared, tatami::make_DelayedBooleanOrVectorHelper<margin_>(std::move(view)));
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
            output->ptr = tatami::make_DelayedUnaryIsometricOp(shared, tatami::make_DelayedBooleanAndScalarHelper(val[0]));
        } else if (op == "|" ){
            output->ptr = tatami::make_DelayedUnaryIsometricOp(shared, tatami::make_DelayedBooleanOrScalarHelper(val[0]));
        } else {
            throw std::runtime_error("unknown delayed boolean operation '" + op + "'");
        }

    } else {
        protectorate[1] = val;
        if (row) {
            set_delayed_boolean_vector<0>(shared, val, op, output->ptr);
        } else {
            set_delayed_boolean_vector<1>(shared, val, op, output->ptr);
        }
    }

    output->original = protectorate; // propagate protection for all child objects by copying references.
    return output;
}

//[[Rcpp::export(rng=false)]]
SEXP apply_delayed_boolean_not(SEXP raw_input) {
    Rtatami::BoundNumericPointer input(raw_input);
    auto output = Rtatami::new_BoundNumericMatrix();
    output->ptr = tatami::make_DelayedUnaryIsometricOp(input->ptr, tatami::DelayedBooleanNotHelper());
    output->original = input->original; // copying the reference to propagate GC protection.
    return output;
}

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
    output->ptr = tatami::make_DelayedUnaryIsometricOp(input->ptr, tatami::DelayedRoundHelper());
    output->original = input->original; // copying the reference to propagate GC protection.
    return output;
}
