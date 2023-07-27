#include "Rtatami.h"
#include "Rcpp.h"
#include "tatami/tatami.hpp"

//[[Rcpp::export(rng=false)]]
SEXP apply_delayed_subset(SEXP raw_input, Rcpp::IntegerVector subset, bool row) {
    Rtatami::BoundNumericPointer input(raw_input);
    auto output = Rtatami::new_BoundNumericMatrix();
    const auto& shared = input->ptr;
    output->original = input->original; // copying the reference to propagate GC protection.

    std::vector<int> resub(subset.begin(), subset.end());
    for (auto& x : resub) {
        --x; 
    } 

    if (row) {
        output->ptr = tatami::make_DelayedSubset<0>(shared, std::move(resub));
    } else {
        output->ptr = tatami::make_DelayedSubset<1>(shared, std::move(resub));
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
