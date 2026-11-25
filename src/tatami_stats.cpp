#include "Rtatami.h"
#include "Rcpp.h"

#include "tatami/tatami.hpp"
#include "tatami_stats/tatami_stats.hpp"
#include "sanisizer/sanisizer.hpp"

#include <vector>
#include <cstddef>
#include <stdexcept>

#include "utils.h"

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector tatami_sums(SEXP raw_input, bool row, int threads) {
    tatami_stats::SumOptions opt;
    opt.num_threads = threads;
    if (threads < 1) {
        throw std::runtime_error("'threads' should be a positive integer");
    }

    Rtatami::BoundNumericPointer input(raw_input);
    const auto NR = input->ptr->nrow();
    const auto NC = input->ptr->ncol();

    auto output = tatami::create_container_of_Index_size<Rcpp::NumericVector>(row ? NR : NC);
    tatami_stats::sum(row, *(input->ptr), static_cast<double*>(output.begin()), opt);
    return output;
}

//[[Rcpp::export(rng=false)]]
Rcpp::List tatami_variances(SEXP raw_input, bool row, int threads) {
    tatami_stats::VarianceOptions opt;
    opt.num_threads = threads;
    if (threads < 1) {
        throw std::runtime_error("'threads' should be a positive integer");
    }

    Rtatami::BoundNumericPointer input(raw_input);
    const auto NR = input->ptr->nrow();
    const auto NC = input->ptr->ncol();

    const auto outlen = (row ? NR : NC);
    auto output_means = tatami::create_container_of_Index_size<Rcpp::NumericVector>(row ? NR : NC);
    auto output_variances = tatami::create_container_of_Index_size<Rcpp::NumericVector>(row ? NR : NC);

    tatami_stats::VarianceBuffers<double> buffers;
    buffers.mean = output_means.begin();
    buffers.variance = output_variances.begin();
    tatami_stats::variance(row, *(input->ptr), buffers, opt);

    return Rcpp::List::create(
        Rcpp::Named("mean") = output_means,
        Rcpp::Named("variance") = output_variances
    );
}

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector tatami_medians(SEXP raw_input, bool row, int threads) {
    tatami_stats::MedianOptions opt;
    opt.num_threads = threads;
    if (threads < 1) {
        throw std::runtime_error("'threads' should be a positive integer");
    }

    Rtatami::BoundNumericPointer input(raw_input);
    const auto NR = input->ptr->nrow();
    const auto NC = input->ptr->ncol();

    auto output = tatami::create_container_of_Index_size<Rcpp::NumericVector>(row ? NR : NC);
    tatami_stats::median(row, *(input->ptr), static_cast<double*>(output.begin()), opt);
    return output;
}

//[[Rcpp::export(rng=false)]]
Rcpp::List tatami_ranges(SEXP raw_input, bool row, int threads) {
    tatami_stats::RangeOptions opt;
    opt.num_threads = threads;
    if (threads < 1) {
        throw std::runtime_error("'threads' should be a positive integer");
    }

    Rtatami::BoundNumericPointer input(raw_input);
    const auto NR = input->ptr->nrow();
    const auto NC = input->ptr->ncol();

    const auto expected = (row ? NR : NC);
    auto minout = tatami::create_container_of_Index_size<Rcpp::NumericVector>(expected);
    auto maxout = tatami::create_container_of_Index_size<Rcpp::NumericVector>(expected);

    tatami_stats::RangeBuffers<double> buf;
    buf.minimum = minout.begin();
    buf.maximum = maxout.begin();

    tatami_stats::range(row, *(input->ptr), buf, opt);
    return Rcpp::List::create(
        Rcpp::Named("minimum") = minout,
        Rcpp::Named("maximum") = maxout
    );
}

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector tatami_nan_counts(SEXP raw_input, bool row, int threads) {
    tatami_stats::CountOptions opt;
    opt.num_threads = threads;
    if (threads < 1) {
        throw std::runtime_error("'threads' should be a positive integer");
    }

    Rtatami::BoundNumericPointer input(raw_input);
    const auto NR = input->ptr->nrow();
    const auto NC = input->ptr->ncol();

    auto output = tatami::create_container_of_Index_size<Rcpp::NumericVector>(row ? NR : NC);
    tatami_stats::count(
        row,
        *(input->ptr),
        static_cast<double*>(output.begin()),
        [](const double val) -> bool { return std::isnan(val); },
        opt
    );
    return output;
}

static std::vector<int> extract_groups(const Rcpp::IntegerVector& group, int num_groups) {
    std::vector<int> group_m1(group.begin(), group.end());
    for (auto& x : group_m1) {
        if (x <= 0) {
            throw std::runtime_error("'group' should contain positive group identifiers");
        } else if (x > num_groups) {
            throw std::runtime_error("entries of 'group' should be less than 'num.groups'");
        }
        --x;
    }
    return group_m1;
}

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector tatami_sums_by_group(SEXP raw_input, Rcpp::IntegerVector group, int num_groups, bool row, int threads) {
    tatami_stats::GroupSumOptions opt;
    opt.num_threads = threads;
    if (threads < 1) {
        throw std::runtime_error("'threads' should be a positive integer");
    }

    Rtatami::BoundNumericPointer input(raw_input);
    const auto NR = input->ptr->nrow();
    const auto NC = input->ptr->ncol();

    auto group_m1 = extract_groups(group, num_groups);
    const auto expected = (row ? NC : NR);
    if (!sanisizer::is_equal(expected, group_m1.size())) {
        throw std::runtime_error("'group' should have length equal to the appropriate dimension extent");
    }

    const auto stride = (row ? NR : NC);
    auto output = create_output_matrix(stride, num_groups);
    auto ptrs = sanisizer::create<std::vector<double*> >(num_groups);
    sanisizer::product<std::size_t>(stride, num_groups); // check that product below is safe.
    for (I<decltype(num_groups)> i = 0; i < num_groups; ++i) {
        ptrs[i] = output.begin() + sanisizer::product_unsafe<std::size_t>(i, stride);
    }

    tatami_stats::group_sum(row, *(input->ptr), group_m1.data(), num_groups, ptrs, opt);

    if (row) {
        return output;
    } else {
        return Rcpp::transpose(output);
    }
}

//[[Rcpp::export(rng=false)]]
Rcpp::List tatami_variances_by_group(SEXP raw_input, Rcpp::IntegerVector group, int num_groups, bool row, int threads) {
    tatami_stats::GroupVarianceOptions opt;
    opt.num_threads = threads;
    if (threads < 1) {
        throw std::runtime_error("'threads' should be a positive integer");
    }

    Rtatami::BoundNumericPointer input(raw_input);
    const auto NR = input->ptr->nrow();
    const auto NC = input->ptr->ncol();

    auto group_m1 = extract_groups(group, num_groups);
    const auto expected = (row ? NC : NR);
    if (!sanisizer::is_equal(expected, group_m1.size())) {
        throw std::runtime_error("'group' should have length equal to the appropriate dimension extent");
    }

    const auto stride = (row ? NR : NC);
    auto moutput = create_output_matrix(stride, num_groups);
    auto voutput = create_output_matrix(stride, num_groups);
    tatami_stats::GroupVarianceBuffers<double> buf;
    buf.mean.reserve(num_groups);
    buf.variance.reserve(num_groups);
    sanisizer::product<std::size_t>(stride, num_groups); // check that all offset calculations are safe.
    for (auto i = static_cast<decltype(num_groups)>(0); i < num_groups; ++i) {
        const auto offset = sanisizer::product_unsafe<std::size_t>(i, stride);
        buf.mean.push_back(moutput.begin() + offset);
        buf.variance.push_back(voutput.begin() + offset);
    }

    tatami_stats::group_variance(row, *(input->ptr), group_m1.data(), num_groups, buf, opt);

    if (row) {
        return Rcpp::List::create(
            Rcpp::Named("mean") = moutput,
            Rcpp::Named("variance") = voutput 
        );
    } else {
        return Rcpp::List::create(
            Rcpp::Named("mean") = Rcpp::transpose(moutput),
            Rcpp::Named("variance") = Rcpp::transpose(voutput)
        );
    }
}

//[[Rcpp::export(rng=false)]]
Rcpp::NumericVector tatami_medians_by_group(SEXP raw_input, Rcpp::IntegerVector group, int num_groups, bool row, int threads) {
    tatami_stats::GroupMedianOptions opt;
    opt.num_threads = threads;
    if (threads < 1) {
        throw std::runtime_error("'threads' should be a positive integer");
    }

    Rtatami::BoundNumericPointer input(raw_input);
    const auto NR = input->ptr->nrow();
    const auto NC = input->ptr->ncol();

    auto group_m1 = extract_groups(group, num_groups);
    const auto expected = (row ? NC : NR);
    if (!sanisizer::is_equal(expected, group_m1.size())) {
        throw std::runtime_error("'group' should have length equal to the appropriate dimension extent");
    }

    const auto stride = (row ? NR : NC);
    auto output = create_output_matrix(stride, num_groups);
    auto ptrs = sanisizer::create<std::vector<double*> >(num_groups);
    sanisizer::product<std::size_t>(stride, num_groups); // check that product inside the loop is safe.
    for (auto i = static_cast<decltype(num_groups)>(0); i < num_groups; ++i) {
        ptrs[i] = output.begin() + sanisizer::product_unsafe<std::size_t>(i, stride);
    }

    tatami_stats::group_median(row, *(input->ptr), group_m1.data(), num_groups, ptrs, opt);

    if (row) {
        return output;
    } else {
        return Rcpp::transpose(output);
    }
}
