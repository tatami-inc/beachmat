#include "Rtatami.h"
#include "Rcpp.h"
#include "tatami/tatami.hpp"

//[[Rcpp::export(rng=false)]]
SEXP apply_delayed_subset(SEXP raw_input, Rcpp::IntegerVector subset, bool row) {
    Rtatami::BoundNumericPointer input(raw_input);
    auto output = Rtatami::new_BoundNumericMatrix();
    const auto& shared = input->ptr;
    output->original = input->original; // copying the reference to propagate GC protection.

    // Is this a contiguous block?
    bool consecutive = true;
    for (size_t i = 1, end = subset.size(); i < end; ++i) {
        if (subset[i] - subset[i - 1] != 1) {
            consecutive = false;
            break;
        }
    }

    if (consecutive) {
        int start = (subset.size() ? subset[0] - 1 : 0);
        int end = (subset.size() ? subset[subset.size() - 1] : 0);
        if (row) {
            output->ptr = tatami::make_DelayedSubsetBlock<0>(shared, start, end);
        } else {
            output->ptr = tatami::make_DelayedSubsetBlock<1>(shared, start, end);
        }
    } else {
        // Otherwise, we get to 0-based indices. This eliminates the need
        // to store 'subset', as we're making our own copy anyway.
        std::vector<int> resub(subset.begin(), subset.end());
        for (auto& x : resub) {
            --x; 
        } 

        if (row) {
            output->ptr = tatami::make_DelayedSubset<0>(shared, std::move(resub));
        } else {
            output->ptr = tatami::make_DelayedSubset<1>(shared, std::move(resub));
        }
    }

    return output;
}

//[[Rcpp::export(rng=false)]]
SEXP apply_delayed_transpose(SEXP raw_input) {
    Rtatami::BoundNumericPointer input(raw_input);
    auto output = Rtatami::new_BoundNumericMatrix();
    output->ptr = tatami::make_DelayedTranspose(input->ptr);
    output->original = input->original; // copying the reference to propagate GC protection.
    return output;
}

//[[Rcpp::export(rng=false)]]
SEXP apply_delayed_bind(Rcpp::List input, bool row) {
    std::vector<std::shared_ptr<tatami::NumericMatrix> > collected;
    collected.reserve(input.size());
    Rcpp::List protectorate(input.size());

    for (size_t i = 0, end = input.size(); i < end; ++i) {
        Rcpp::RObject current = input[i];
        Rtatami::BoundNumericPointer curptr(current);
        protectorate[i] = curptr->original;
        collected.push_back(curptr->ptr);
    }

    auto output = Rtatami::new_BoundNumericMatrix();
    if (row) {
        output->ptr = tatami::make_DelayedBind<0>(std::move(collected));
    } else {
        output->ptr = tatami::make_DelayedBind<1>(std::move(collected));
    }

    output->original = protectorate; // propagate protection for all child objects by copying references.
    return output;
}

//[[Rcpp::export(rng=false)]]
SEXP apply_delayed_addition(SEXP raw_input, Rcpp::NumericVector val, bool row) {
    Rtatami::BoundNumericPointer input(raw_input);
    const auto& shared = input->ptr;
    Rcpp::List protectorate(2);
    protectorate[0] = input->original;

    auto output = Rtatami::new_BoundNumericMatrix();
    if (val.size() == 1) {
        output->ptr = tatami::make_DelayedIsometricOp(shared, tatami::DelayedAddScalarHelper(val[0]));
    } else {
        protectorate[1] = val;
        tatami::ArrayView<double> view(static_cast<const double*>(val.begin()), val.size());
        if (row) {
            output->ptr = tatami::make_DelayedIsometricOp(shared, tatami::make_DelayedAddVectorHelper<0>(std::move(view)));
        } else {
            output->ptr = tatami::make_DelayedIsometricOp(shared, tatami::make_DelayedAddVectorHelper<1>(std::move(view)));
        }
    }

    output->original = protectorate; // propagate protection for all child objects by copying references.
    return output;
}

//[[Rcpp::export(rng=false)]]
SEXP apply_delayed_multiplication(SEXP raw_input, Rcpp::NumericVector val, bool row) {
    Rtatami::BoundNumericPointer input(raw_input);
    const auto& shared = input->ptr;
    Rcpp::List protectorate(2);
    protectorate[0] = input->original;

    auto output = Rtatami::new_BoundNumericMatrix();
    if (val.size() == 1) {
        output->ptr = tatami::make_DelayedIsometricOp(shared, tatami::DelayedMultiplyScalarHelper(val[0]));
    } else {
        protectorate[1] = val;
        tatami::ArrayView<double> view(static_cast<const double*>(val.begin()), val.size());
        if (row) {
            output->ptr = tatami::make_DelayedIsometricOp(shared, tatami::make_DelayedMultiplyVectorHelper<0>(std::move(view)));
        } else {
            output->ptr = tatami::make_DelayedIsometricOp(shared, tatami::make_DelayedMultiplyVectorHelper<1>(std::move(view)));
        }
    }

    output->original = protectorate; // propagate protection for all child objects by copying references.
    return output;
}

//[[Rcpp::export(rng=false)]]
SEXP apply_delayed_subtraction(SEXP raw_input, Rcpp::NumericVector val, bool right, bool row) {
    Rtatami::BoundNumericPointer input(raw_input);
    const auto& shared = input->ptr;
    Rcpp::List protectorate(2);
    protectorate[0] = input->original;

    auto output = Rtatami::new_BoundNumericMatrix();
    if (val.size() == 1) {
        if (right) {
            output->ptr = tatami::make_DelayedIsometricOp(shared, tatami::DelayedSubtractScalarHelper<true>(val[0]));
        } else {
            output->ptr = tatami::make_DelayedIsometricOp(shared, tatami::DelayedSubtractScalarHelper<false>(val[0]));
        }

    } else {
        protectorate[1] = val;        
        tatami::ArrayView<double> view(static_cast<const double*>(val.begin()), val.size());
        if (right) {
            if (row) {
                output->ptr = tatami::make_DelayedIsometricOp(shared, tatami::make_DelayedSubtractVectorHelper<true, 0>(std::move(view)));
            } else {
                output->ptr = tatami::make_DelayedIsometricOp(shared, tatami::make_DelayedSubtractVectorHelper<true, 1>(std::move(view)));
            }
        } else {
            if (row) {
                output->ptr = tatami::make_DelayedIsometricOp(shared, tatami::make_DelayedSubtractVectorHelper<false, 0>(std::move(view)));
            } else {
                output->ptr = tatami::make_DelayedIsometricOp(shared, tatami::make_DelayedSubtractVectorHelper<false, 1>(std::move(view)));
            }
        }
    }

    output->original = protectorate; // propagate protection for all child objects by copying references.
    return output;
}

//[[Rcpp::export(rng=false)]]
SEXP apply_delayed_division(SEXP raw_input, Rcpp::NumericVector val, bool right, bool row) {
    Rtatami::BoundNumericPointer input(raw_input);
    const auto& shared = input->ptr;
    Rcpp::List protectorate(2);
    protectorate[0] = input->original;

    auto output = Rtatami::new_BoundNumericMatrix();
    if (val.size() == 1) {
        if (right) {
            output->ptr = tatami::make_DelayedIsometricOp(shared, tatami::DelayedDivideScalarHelper<true>(val[0]));
        } else {
            output->ptr = tatami::make_DelayedIsometricOp(shared, tatami::DelayedDivideScalarHelper<false>(val[0]));
        }

    } else {
        protectorate[1] = val;        
        tatami::ArrayView<double> view(static_cast<const double*>(val.begin()), val.size());
        if (right) {
            if (row) {
                output->ptr = tatami::make_DelayedIsometricOp(shared, tatami::make_DelayedDivideVectorHelper<true, 0>(std::move(view)));
            } else {
                output->ptr = tatami::make_DelayedIsometricOp(shared, tatami::make_DelayedDivideVectorHelper<true, 1>(std::move(view)));
            }
        } else {
            if (row) {
                output->ptr = tatami::make_DelayedIsometricOp(shared, tatami::make_DelayedDivideVectorHelper<false, 0>(std::move(view)));
            } else {
                output->ptr = tatami::make_DelayedIsometricOp(shared, tatami::make_DelayedDivideVectorHelper<false, 1>(std::move(view)));
            }
        }
    }

    output->original = protectorate; // propagate protection for all child objects by copying references.
    return output;
}

//[[Rcpp::export(rng=false)]]
SEXP apply_delayed_log(SEXP raw_input, double base) {
    Rtatami::BoundNumericPointer input(raw_input);
    auto output = Rtatami::new_BoundNumericMatrix();
    output->ptr = tatami::make_DelayedIsometricOp(input->ptr, tatami::DelayedLogHelper(base));
    output->original = input->original; // copying the reference to propagate GC protection.
    return output;
}

//[[Rcpp::export(rng=false)]]
SEXP apply_delayed_log1p(SEXP raw_input) {
    Rtatami::BoundNumericPointer input(raw_input);
    auto output = Rtatami::new_BoundNumericMatrix();
    output->ptr = tatami::make_DelayedIsometricOp(input->ptr, tatami::DelayedLog1pHelper<>());
    output->original = input->original; // copying the reference to propagate GC protection.
    return output;
}

//[[Rcpp::export(rng=false)]]
SEXP apply_delayed_abs(SEXP raw_input) {
    Rtatami::BoundNumericPointer input(raw_input);
    auto output = Rtatami::new_BoundNumericMatrix();
    output->ptr = tatami::make_DelayedIsometricOp(input->ptr, tatami::DelayedAbsHelper<>());
    output->original = input->original; // copying the reference to propagate GC protection.
    return output;
}

//[[Rcpp::export(rng=false)]]
SEXP apply_delayed_sqrt(SEXP raw_input) {
    Rtatami::BoundNumericPointer input(raw_input);
    auto output = Rtatami::new_BoundNumericMatrix();
    output->ptr = tatami::make_DelayedIsometricOp(input->ptr, tatami::DelayedSqrtHelper<>());
    output->original = input->original; // copying the reference to propagate GC protection.
    return output;
}

//[[Rcpp::export(rng=false)]]
SEXP apply_delayed_round(SEXP raw_input) {
    Rtatami::BoundNumericPointer input(raw_input);
    auto output = Rtatami::new_BoundNumericMatrix();
    output->ptr = tatami::make_DelayedIsometricOp(input->ptr, tatami::DelayedRoundHelper<>());
    output->original = input->original; // copying the reference to propagate GC protection.
    return output;
}

//[[Rcpp::export(rng=false)]]
SEXP apply_delayed_exp(SEXP raw_input) {
    Rtatami::BoundNumericPointer input(raw_input);
    auto output = Rtatami::new_BoundNumericMatrix();
    output->ptr = tatami::make_DelayedIsometricOp(input->ptr, tatami::DelayedExpHelper<>());
    output->original = input->original; // copying the reference to propagate GC protection.
    return output;
}
