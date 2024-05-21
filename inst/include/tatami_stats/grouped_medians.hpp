#ifndef TATAMI_STATS_GROUPED_MEDIANS_HPP
#define TATAMI_STATS_GROUPED_MEDIANS_HPP

#include "utils.hpp"
#include "tatami/tatami.hpp"
#include "medians.hpp"
#include <vector>
#include <algorithm>

/**
 * @file grouped_medians.hpp
 *
 * @brief Compute group-wise medians from a `tatami::Matrix`.
 */

namespace tatami_stats {

/**
 * @brief Functions for computing dimension-wise grouped medians.
 * @namespace tatami_stats::grouped_medians
 */
namespace grouped_medians {

/**
 * @brief Grouped median calculation options.
 */
struct Options {
    /**
     * Whether to check for NaNs in the input, and skip them.
     * If false, NaNs are assumed to be absent, and the behavior of the median calculation in the presence of NaNs is undefined.
     */
    bool skip_nan = false;

    /**
     * Number of threads to use when computing medians across a `tatami::Matrix`.
     */
    int num_threads = 1;
};

/**
 * Compute per-group medians for each element of a chosen dimension of a `tatami::Matrix`.
 *
 * @tparam Value_ Type of the matrix value, should be numeric.
 * @tparam Index_ Type of the row/column indices.
 * @tparam Group_ Type of the group assignments for each column.
 * @tparam GroupSizes_ Vector-like class that has `size()` and `[` methods.
 * @tparam Output_ Type of the output value.
 * This should be floating-point to store potential averages.
 *
 * @param row Whether to compute medians for the rows.
 * @param p Pointer to a `tatami::Matrix`.
 * @param[in] group Pointer to an array of length equal to the number of columns (if `row = true`) or rows (otherwise).
 * Each value should be an integer that specifies the group assignment.
 * Values should lie in `[0, N)` where `N` is the number of unique groups.
 * @param group_sizes Vector-like object of length `N`, specifying the number of columns assigned to each group.
 * This can be created by calling `tatami_stats::tabulate_groups()` on `group`.
 * @param[out] output Pointer to an array of pointers of length equal to the number of groups.
 * Each inner pointer should reference an array of length equal to the number of rows (if `row = true`) or columns (otherwise).
 * On output, this will contain the row/column medians for each group (indexed according to the assignment in `group`).
 * @param mopt Median calculation options.
 */
template<typename Value_, typename Index_, typename Group_, class GroupSizes_, typename Output_>
void apply(bool row, const tatami::Matrix<Value_, Index_>* p, const Group_* group, const GroupSizes_& group_sizes, Output_** output, const Options& mopt) {
    Index_ dim = (row ? p->nrow() : p->ncol());
    Index_ otherdim = (row ? p->ncol() : p->nrow());

    tatami::parallelize([&](int, Index_ start, Index_ len) -> void {
        std::vector<Value_> xbuffer(otherdim);

        size_t ngroups = group_sizes.size();
        std::vector<std::vector<double> > workspace(ngroups);
        for (size_t g = 0; g < ngroups; ++g) {
            workspace[g].reserve(group_sizes[g]);
        }

        if (p->sparse()) {
            tatami::Options opt;
            opt.sparse_ordered_index = false;

            auto ext = tatami::consecutive_extractor<true>(p, row, start, len, opt);
            std::vector<Index_> ibuffer(otherdim);
            for (Index_ i = 0; i < len; ++i) {
                auto range = ext->fetch(xbuffer.data(), ibuffer.data());
                for (Index_ j = 0; j < range.number; ++j) {
                    workspace[group[range.index[j]]].push_back(range.value[j]);
                }

                for (size_t g = 0; g < ngroups; ++g) {
                    auto& w = workspace[g];
                    output[g][i + start] = medians::direct(w.data(), w.size(), static_cast<size_t>(group_sizes[g]), mopt.skip_nan);
                    w.clear();
                }
            }

        } else {
            auto ext = tatami::consecutive_extractor<false>(p, row, start, len);
            for (Index_ i = 0; i < len; ++i) {
                auto ptr = ext->fetch(xbuffer.data());
                for (Index_ j = 0; j < otherdim; ++j) {
                    workspace[group[j]].push_back(ptr[j]);
                }

                for (size_t g = 0; g < ngroups; ++g) {
                    auto& w = workspace[g];
                    output[g][i + start] = medians::direct(w.data(), w.size(), mopt.skip_nan);
                    w.clear();
                }
            }
        }
    }, dim, mopt.num_threads);
}

/**
 * Wrapper around `apply()` for row-wise grouped medians.
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
 * @param mopt Median calculation options.
 *
 * @return Vector of length equal to the number of groups.
 * Each entry is a vector of length equal to the number of rows, containing the row-wise medians for the corresponding group.
 */
template<typename Output_ = double, typename Value_, typename Index_, typename Group_>
std::vector<std::vector<Output_> > by_row(const tatami::Matrix<Value_, Index_>* p, const Group_* group, const Options& mopt) {
    size_t mydim = p->nrow();
    auto group_sizes = tabulate_groups(group, p->ncol());

    std::vector<std::vector<Output_> > output(group_sizes.size());
    std::vector<Output_*> ptrs;
    ptrs.reserve(output.size());
    for (auto& o : output) {
        o.resize(mydim);
        ptrs.push_back(o.data());
    }

    apply(true, p, group, group_sizes, ptrs.data(), mopt);
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
 * Each entry is a vector of length equal to the number of rows, containing the row-wise medians for the corresponding group.
 */
template<typename Output_ = double, typename Value_, typename Index_, typename Group_>
std::vector<std::vector<Output_> > by_row(const tatami::Matrix<Value_, Index_>* p, const Group_* group) {
    return by_row(p, group, Options());
}

/**
 * Wrapper around `apply()` for column-wise grouped medians.
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
 * @param mopt Median calculation options.
 *
 * @return Vector of length equal to the number of groups.
 * Each entry is a vector of length equal to the number of columns, containing the column-wise medians for the corresponding group.
 */
template<typename Output_ = double, typename Value_, typename Index_, typename Group_>
std::vector<std::vector<Output_> > by_column(const tatami::Matrix<Value_, Index_>* p, const Group_* group, const Options& mopt) {
    size_t mydim = p->ncol();
    auto group_sizes = tabulate_groups(group, p->nrow());

    std::vector<std::vector<Output_> > output(group_sizes.size());
    std::vector<Output_*> ptrs;
    ptrs.reserve(output.size());
    for (auto& o : output) {
        o.resize(mydim);
        ptrs.push_back(o.data());
    }

    apply(false, p, group, group_sizes, ptrs.data(), mopt);
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
 * Each entry is a vector of length equal to the number of columns, containing the column-wise medians for the corresponding group.
 */
template<typename Output_ = double, typename Value_, typename Index_, typename Group_>
std::vector<std::vector<Output_> > by_column(const tatami::Matrix<Value_, Index_>* p, const Group_* group) {
    return by_column(p, group, Options());
}

}

}

#endif
