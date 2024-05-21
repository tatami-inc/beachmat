#ifndef TATAMI_STATS_GROUPED_SUMS_HPP
#define TATAMI_STATS_GROUPED_SUMS_HPP

#include "utils.hpp"
#include "tatami/tatami.hpp"
#include "sums.hpp"
#include <vector>
#include <algorithm>

/**
 * @file grouped_sums.hpp
 *
 * @brief Compute group-wise sums from a `tatami::Matrix`.
 */

namespace tatami_stats {

/**
 * @brief Functions for computing dimension-wise grouped sums.
 * @namespace tatami_stats::grouped_sums
 */
namespace grouped_sums {

/**
 * @brief Grouped summation options.
 */
struct Options {
    /**
     * Whether to check for NaNs in the input, and skip them.
     * If false, NaNs are assumed to be absent, and the behavior of the summation in the presence of NaNs is undefined.
     */
    bool skip_nan = false;

    /**
     * Number of threads to use when computing sums across a `tatami::Matrix`.
     */
    int num_threads = 1;
};

/**
 * Compute per-group sums for each element of a chosen dimension of a `tatami::Matrix`.
 *
 * @tparam Value_ Type of the matrix value, should be numeric.
 * @tparam Index_ Type of the row/column indices.
 * @tparam Group_ Type of the group assignments for each column.
 * @tparam Output_ Type of the output value.
 * This should be floating-point to store potential averages.
 *
 * @param row Whether to compute sums for the rows.
 * @param p Pointer to a `tatami::Matrix`.
 * @param[in] group Pointer to an array of length equal to the number of columns (if `row = true`) or rows (otherwise).
 * Each value should be an integer that specifies the group assignment.
 * Values should lie in `[0, N)` where `N` is the number of unique groups.
 * @param num_groups Number of groups, i.e., `N`.
 * This can be determined by calling `tatami_stats::total_groups()` on `group`.
 * @param[out] output Pointer to an array of pointers of length equal to the number of groups.
 * Each inner pointer should reference an array of length equal to the number of rows (if `row = true`) or columns (otherwise).
 * On output, this will contain the row/column sums for each group (indexed according to the assignment in `group`).
 * @param sopt Summation options.
 */
template<typename Value_, typename Index_, typename Group_, typename Output_>
void apply(bool row, const tatami::Matrix<Value_, Index_>* p, const Group_* group, size_t num_groups, Output_** output, const Options& sopt) {
    Index_ dim = (row ? p->nrow() : p->ncol());
    Index_ otherdim = (row ? p->ncol() : p->nrow());

    if (p->sparse()) {
        if (p->prefer_rows() == row) {
            tatami::parallelize([&](int, Index_ start, Index_ len) -> void {
                auto ext = tatami::consecutive_extractor<true>(p, row, start, len);
                std::vector<Value_> xbuffer(otherdim);
                std::vector<Index_> ibuffer(otherdim);
                std::vector<Output_> tmp(num_groups);

                for (Index_ i = 0; i < len; ++i) {
                    auto range = ext->fetch(xbuffer.data(), ibuffer.data());
                    std::fill(tmp.begin(), tmp.end(), static_cast<Output_>(0));

                    if (sopt.skip_nan) {
                        for (int j = 0; j < range.number; ++j) {
                            auto val = range.value[j];
                            if (!std::isnan(val)) {
                                tmp[group[range.index[j]]] += val;
                            }
                        }
                    } else {
                        for (int j = 0; j < range.number; ++j) {
                            tmp[group[range.index[j]]] += range.value[j];
                        }
                    }

                    for (size_t g = 0; g < num_groups; ++g) {
                        output[g][i + start] = tmp[g];
                    }
                }
            }, dim, sopt.num_threads);

        } else {
            // Order within each observed vector doesn't affect numerical
            // precision of the outcome, as addition order for each objective
            // vector is already well-defined for a running calculation.
            tatami::Options opt;
            opt.sparse_ordered_index = false; 

            tatami::parallelize([&](size_t thread, Index_ start, Index_ len) -> void {
                std::vector<sums::RunningSparse<Output_, Value_, Index_> > runners;
                runners.reserve(num_groups);
                std::vector<LocalOutputBuffer<Output_> > local_output;
                local_output.reserve(num_groups);

                for (size_t g = 0; g < num_groups; ++g) {
                    local_output.emplace_back(thread, start, len, output[g]);
                    runners.emplace_back(local_output.back().data(), sopt.skip_nan, start);
                }

                auto ext = tatami::consecutive_extractor<true>(p, !row, 0, otherdim, start, len, opt);
                std::vector<Value_> xbuffer(len);
                std::vector<Index_> ibuffer(len);

                for (int i = 0; i < otherdim; ++i) {
                    auto range = ext->fetch(xbuffer.data(), ibuffer.data());
                    runners[group[i]].add(range.value, range.index, range.number);
                }

                for (size_t g = 0; g < num_groups; ++g) {
                    local_output[g].transfer();
                }
            }, dim, sopt.num_threads);
        }

    } else {
        if (p->prefer_rows() == row) {
            tatami::parallelize([&](int, Index_ start, Index_ len) -> void {
                auto ext = tatami::consecutive_extractor<false>(p, row, start, len);
                std::vector<Value_> xbuffer(otherdim);
                std::vector<Output_> tmp(num_groups);

                for (Index_ i = 0; i < len; ++i) {
                    auto ptr = ext->fetch(xbuffer.data());
                    std::fill(tmp.begin(), tmp.end(), static_cast<Output_>(0));

                    if (sopt.skip_nan) {
                        for (Index_ j = 0; j < otherdim; ++j) {
                            auto val = ptr[j];
                            if (!std::isnan(val)) {
                                tmp[group[j]] += val;
                            }
                        }
                    } else {
                        for (Index_ j = 0; j < otherdim; ++j) {
                            tmp[group[j]] += ptr[j];
                        }
                    }

                    for (size_t g = 0; g < num_groups; ++g) {
                        output[g][i + start] = tmp[g];
                    }
                }
            }, dim, sopt.num_threads);

        } else {
            tatami::parallelize([&](size_t thread, Index_ start, Index_ len) -> void {
                std::vector<sums::RunningDense<Output_, Value_, Index_> > runners;
                runners.reserve(num_groups);
                std::vector<LocalOutputBuffer<Output_> > local_output;
                local_output.reserve(num_groups);

                for (size_t g = 0; g < num_groups; ++g) {
                    local_output.emplace_back(thread, start, len, output[g]);
                    runners.emplace_back(len, local_output.back().data(), sopt.skip_nan);
                }

                std::vector<double> xbuffer(len);
                auto ext = tatami::consecutive_extractor<false>(p, !row, 0, otherdim, start, len);

                for (int i = 0; i < otherdim; ++i) {
                    auto ptr = ext->fetch(xbuffer.data());
                    runners[group[i]].add(ptr);
                }

                for (size_t g = 0; g < num_groups; ++g) {
                    local_output[g].transfer();
                }
            }, dim, sopt.num_threads);
        }
    }
}

/**
 * Wrapper around `apply()` for row-wise grouped sums.
 *
 * @tparam Output_ Type of the output.
 * @tparam Value_ Type of the matrix value.
 * @tparam Index_ Type of the row/column indices.
 * @tparam Group_ Type of the group assignments for each row.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param[in] group Pointer to an array of length equal to the number of columns.
 * Each value should be an integer that specifies the group assignment.
 * Values should lie in `[0, N)` where `N` is the number of unique groups.
 * @param sopt Summation options.
 *
 * @return Vector of length equal to the number of groups.
 * Each entry is a vector of length equal to the number of rows, containing the row-wise sums for the corresponding group.
 */
template<typename Output_ = double, typename Value_, typename Index_, typename Group_>
std::vector<std::vector<Output_> > by_row(const tatami::Matrix<Value_, Index_>* p, const Group_* group, const Options& sopt) {
    size_t mydim = p->nrow();
    auto ngroup = total_groups(group, p->ncol());

    std::vector<std::vector<Output_> > output(ngroup);
    std::vector<Output_*> ptrs;
    ptrs.reserve(output.size());
    for (auto& o : output) {
        o.resize(mydim);
        ptrs.push_back(o.data());
    }

    apply(true, p, group, ngroup, ptrs.data(), sopt);
    return output;
}

/**
 * Overload with default options.
 *
 * @tparam Output_ Type of the output.
 * @tparam Value_ Type of the matrix value.
 * @tparam Index_ Type of the row/column indices.
 * @tparam Group_ Type of the group assignments for each column.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param[in] group Pointer to an array of length equal to the number of columns.
 * Each value should be an integer that specifies the group assignment.
 * Values should lie in `[0, N)` where `N` is the number of unique groups.
 *
 * @return Vector of length equal to the number of groups.
 * Each entry is a vector of length equal to the number of rows, containing the row-wise sums for the corresponding group.
 */
template<typename Output_ = double, typename Value_, typename Index_, typename Group_>
std::vector<std::vector<Output_> > by_row(const tatami::Matrix<Value_, Index_>* p, const Group_* group) {
    return by_row(p, group, Options());
}

/**
 * Wrapper around `apply()` for column-wise grouped sums.
 *
 * @tparam Output_ Type of the output.
 * @tparam Value_ Type of the matrix value.
 * @tparam Index_ Type of the column/column indices.
 * @tparam Group_ Type of the group assignments for each column.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param[in] group Pointer to an array of length equal to the number of rows.
 * Each value should be an integer that specifies the group assignment.
 * Values should lie in `[0, N)` where `N` is the number of unique groups.
 * @param sopt Summation options.
 *
 * @return Vector of length equal to the number of groups.
 * Each entry is a vector of length equal to the number of columns, containing the column-wise sums for the corresponding group.
 */
template<typename Output_ = double, typename Value_, typename Index_, typename Group_>
std::vector<std::vector<Output_> > by_column(const tatami::Matrix<Value_, Index_>* p, const Group_* group, const Options& sopt) {
    size_t mydim = p->ncol();
    auto ngroup = total_groups(group, p->nrow());

    std::vector<std::vector<Output_> > output(ngroup);
    std::vector<Output_*> ptrs;
    ptrs.reserve(output.size());
    for (auto& o : output) {
        o.resize(mydim);
        ptrs.push_back(o.data());
    }

    apply(false, p, group, ngroup, ptrs.data(), sopt);
    return output;
}

/**
 * Overload with default options.
 *
 * @tparam Output_ Type of the output.
 * @tparam Value_ Type of the matrix value.
 * @tparam Index_ Type of the column/column indices.
 * @tparam Group_ Type of the group assignments for each column.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param[in] group Pointer to an array of length equal to the number of rows.
 * Each value should be an integer that specifies the group assignment.
 * Values should lie in `[0, N)` where `N` is the number of unique groups.
 *
 * @return Vector of length equal to the number of groups.
 * Each entry is a vector of length equal to the number of columns, containing the column-wise sums for the corresponding group.
 */
template<typename Output_ = double, typename Value_, typename Index_, typename Group_>
std::vector<std::vector<Output_> > by_column(const tatami::Matrix<Value_, Index_>* p, const Group_* group) {
    return by_column(p, group, Options());
}

}

}

#endif
