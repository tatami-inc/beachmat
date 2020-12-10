#include "Rcpp.h"
#include <vector>
#include <algorithm>

// [[Rcpp::export(rng=false)]]
Rcpp::List fragment_sparse_rows(Rcpp::IntegerVector i, Rcpp::IntegerVector p, Rcpp::IntegerVector limits) {
    size_t ncolsp1 = p.size();
    size_t nchunks = limits.size();

    std::vector<Rcpp::IntegerVector> starti(nchunks);
    std::vector<Rcpp::IntegerVector> newp(nchunks);
    for (size_t l = 0; l < nchunks; ++l) {
        starti[l] = Rcpp::IntegerVector(ncolsp1 - 1);
        newp[l] = Rcpp::IntegerVector(ncolsp1);
    }

    int counter = 0;
    auto iIt = i.begin();

    for (size_t px = 1; px < ncolsp1; ++px) {
        auto colend = p[px];
        for (size_t l = 0; l < nchunks; ++l) {
            auto chunkend = limits[l];

            auto& start = starti[l][px - 1];
            start = counter;
            while (iIt != i.end() && counter < colend && *iIt < chunkend) {
                ++iIt;
                ++counter;
            }

            auto& curp = newp[l];
            curp[px] = curp[px - 1] + counter - start;
        }
    }

    Rcpp::List output(nchunks);
    for (size_t l = 0; l < nchunks; ++l) {
        output[l] = Rcpp::List::create(starti[l], newp[l], R_NilValue, R_NilValue);
    }

    return output;
}

// [[Rcpp::export(rng=false)]]
Rcpp::IntegerVector sparse_subset_index (Rcpp::IntegerVector starts, Rcpp::IntegerVector newp) {
    size_t ncols = starts.size();
    size_t total = newp[ncols];
    Rcpp::IntegerVector subset(total);
    auto sIt = subset.begin();
    for (size_t px = 1; px <= ncols; ++px) {
        auto delta = newp[px] - newp[px-1];
        std::iota(sIt, sIt + delta, starts[px-1] + 1);
        sIt += delta;
    }
    return subset;
}
