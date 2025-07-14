#include "Rtatami.h"
#include "Rcpp.h"
#include "tatami/tatami.hpp"

#include <vector>
#include <memory>

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

    output->ptr = tatami::make_DelayedSubset(shared, std::move(resub), row);
    return output;
}

//[[Rcpp::export(rng=false)]]
SEXP apply_delayed_transpose(SEXP raw_input) {
    Rtatami::BoundNumericPointer input(raw_input);
    auto output = Rtatami::new_BoundNumericMatrix();
    output->ptr.reset(new tatami::DelayedTranspose<double, int>(input->ptr));
    output->original = input->original; // copying the reference to propagate GC protection.
    return output;
}

//[[Rcpp::export(rng=false)]]
SEXP apply_delayed_bind(Rcpp::List input, bool row) {
    std::vector<std::shared_ptr<tatami::NumericMatrix> > collected;
    collected.reserve(input.size());
    Rcpp::List protectorate(input.size());

    for (decltype(input.size()) i = 0, end = input.size(); i < end; ++i) {
        Rcpp::RObject current = input[i];
        Rtatami::BoundNumericPointer curptr(current);
        protectorate[i] = curptr->original;
        collected.push_back(curptr->ptr);
    }

    auto output = Rtatami::new_BoundNumericMatrix();
    output->ptr.reset(new tatami::DelayedBind<double, int>(std::move(collected), row));
    output->original = protectorate; // propagate protection for all child objects by copying references.
    return output;
}
