#ifndef TATAMI_STATS_RANGES_HPP
#define TATAMI_STATS_RANGES_HPP

#include "tatami/tatami.hpp"
#include "utils.hpp"

#include <vector>
#include <algorithm>
#include <type_traits>

/**
 * @file ranges.hpp
 *
 * @brief Compute row and column ranges from a `tatami::Matrix`.
 */

namespace tatami_stats {

/**
 * @brief Functions for computing dimension-wise ranges.
 * @namespace tatami_stats::ranges
 */
namespace ranges {

/**
 * @brief Range calculation options.
 */
struct Options {
    /**
     * Whether to check for NaNs in the input, and skip them.
     * If false, NaNs are assumed to be absent, and the behavior of the range calculation in the presence of NaNs is undefined.
     */
    bool skip_nan = false;

    /**
     * Number of threads to use when computing ranges across a `tatami::Matrix`.
     */
    int num_threads = 1;
};

/**
 * @cond
 */
namespace internal {

template<bool minimum_, typename Value_>
constexpr auto choose_placeholder() {
    if constexpr(minimum_) {
        // Placeholder value 'x' is such that 'x > y' is always true for any non-NaN 'y'.
        if constexpr(std::numeric_limits<Value_>::has_infinity) {
            return std::numeric_limits<Value_>::infinity();
        } else {
            return std::numeric_limits<Value_>::max();
        }
    } else {
        // Placeholder value 'x' is such that 'x < y' is always true for any non-NaN 'y'.
        if constexpr(std::numeric_limits<Value_>::has_infinity) {
            return -std::numeric_limits<Value_>::infinity();
        } else {
            return std::numeric_limits<Value_>::lowest();
        }
    }
}

template<bool minimum_, typename Output_, typename Value_>
bool is_better(Output_ best, Value_ alt) {
    if constexpr(minimum_) {
        return best > static_cast<Output_>(alt);
    } else {
        return best < static_cast<Output_>(alt);
    }
}

}
/**
 * @endcond
 */

/**
 * Directly compute the minimum or maximum of a dense array.
 *
 * @tparam minimum_ Whether to compute the minimum.
 * If false, the maximum is computed instead.
 * @tparam Value_ Type of the input data.
 * @tparam Index_ Type of the row/column indices.
 *
 * @param[in] ptr Pointer to an array of values of length `num`.
 * @param num Size of the array.
 * @param skip_nan See `Options::skip_nan` for details.
 *
 * @return The minimum or maximum value, depending on `minimum_`.
 * If `num = 0` or (if `skip_nan_ = true`) there are no non-NaN values, a placeholder value is returned instead
 * that is never less than (if `minimum_ true`) or greater than (otherwise) any non-NaN value of type `Value_`.
 */
template<bool minimum_, typename Value_, typename Index_>
Value_ direct(const Value_* ptr, Index_ num, bool skip_nan) {
    if (skip_nan) {
        auto current = internal::choose_placeholder<minimum_, Value_>(); 
        for (Index_ i = 0; i < num; ++i) {
            auto val = ptr[i];
            if (internal::is_better<minimum_>(current, val)) { // no need to explicitly handle NaNs, as any comparison with NaNs is always false.
                current = val;
            }
        }
        return current;

    } else if (num) {
        if constexpr(minimum_) {
            return *std::min_element(ptr, ptr + num);
        } else {
            return *std::max_element(ptr, ptr + num);
        }

    } else {
        return internal::choose_placeholder<minimum_, Value_>(); 
    }
}

/**
 * Compute the extremes of a sparse array.
 *
 * @tparam minimum_ Whether to compute the minimum.
 * If false, the maximum is computed instead.
 * @tparam Value_ Type of the input data.
 * @tparam Index_ Type of the row/column indices.
 *
 * @param[in] value Pointer to an array of values of length `num`.
 * @param num_nonzero Length of the array pointed to by `value`.
 * @param num_all Total number of values in the dataset, including the zeros not in `value`.
 * This should be greater than or equal to `num_nonzero`.
 * @param skip_nan See `Options::skip_nan` for details.
 *
 * @return The minimum or maximum value, depending on `minimum_`.
 * If `num_all = 0` or (if `skip_nan_ = true`) there are no non-NaN values, a placeholder value is returned instead
 * that is never less than (if `minimum_ true`) or greater than (otherwise) any non-NaN value of type `Value_`.
 */
template<bool minimum_, typename Value_, typename Index_>
Value_ direct(const Value_* value, Index_ num_nonzero, Index_ num_all, bool skip_nan) {
    if (num_nonzero) {
        auto candidate = direct<minimum_>(value, num_nonzero, skip_nan);
        if (num_nonzero < num_all && internal::is_better<minimum_>(candidate, 0)) {
            candidate = 0;
        }
        return candidate;
    } else if (num_all) {
        return 0;
    } else {
        return internal::choose_placeholder<minimum_, Value_>();
    }
}

/**
 * @brief Running minima/maxima from dense data.
 *
 * This considers a scenario with a set of equilength "objective" vectors [V1, V2, V3, ..., Vn],
 * but data are only available for "observed" vectors [P1, P2, P3, ..., Pm],
 * where Pi[j] contains the i-th element of objective vector Vj.
 * The idea is to repeatedly call `add()` for `ptr` corresponding to observed vectors from 0 to m - 1,
 * which computes the running minimum/maximum for each objective vector at each invocation.
 *
 * @tparam minimum_ Whether to compute the minimum.
 * @tparam Output_ Type of the output data.
 * @tparam Value_ Type of the input data.
 * @tparam Index_ Type of the row/column indices.
 */
template<bool minimum_, typename Output_, typename Value_, typename Index_>
class RunningDense {
public:
    /**
     * @param num Number of objective vectors, i.e., n.
     * @param[out] store Pointer to an output array of length `num`.
     * After `finish()` is called, this will contain the minimum/maximum for each objective vector.
     * @param skip_nan See `Options::skip_nan` for details.
     */
    RunningDense(Index_ num, Output_* store, bool skip_nan) : my_num(num), my_store(store), my_skip_nan(skip_nan) {}

    /**
     * Add the next observed vector to the running min/max calculation.
     * @param[in] ptr Pointer to an array of values of length `num`, corresponding to an observed vector.
     */
    void add(const Value_* ptr) {
        if (my_init) {
            my_init = false;
            if (my_skip_nan) {
                for (Index_ i = 0; i < my_num; ++i, ++ptr) {
                    auto val = *ptr;
                    if (std::isnan(val)) {
                        my_store[i] = internal::choose_placeholder<minimum_, Value_>();
                    } else {
                        my_store[i] = val;
                    }
                }
            } else {
                std::copy_n(ptr, my_num, my_store);
            }

        } else {
            for (Index_ i = 0; i < my_num; ++i, ++ptr) {
                auto val = *ptr;
                if (internal::is_better<minimum_>(my_store[i], val)) { // this should implicitly skip NaNs, any NaN comparison will be false.
                    my_store[i] = val;
                }
            }
        }
    }

    /**
     * Finish the running calculation once all observed vectors have been passed to `add()`. 
     */
    void finish() {
        if (my_init) {
            std::fill_n(my_store, my_num, internal::choose_placeholder<minimum_, Value_>());
        }
    }

private:
    bool my_init = true;
    Index_ my_num;
    Output_* my_store;
    bool my_skip_nan;
};

/**
 * @brief Running minima/maxima from sparse data.
 *
 * Compute running minima and maximuma from sparse data. 
 * This does the same as `RunningDense` but for sparse observed vectors.
 *
 * @tparam minimum_ Whether to compute the minimum.
 * @tparam Output_ Type of the output data.
 * @tparam Value_ Type of the input value.
 * @tparam Index_ Type of the row/column indices.
 */
template<bool minimum_, typename Output_, typename Value_, typename Index_>
class RunningSparse {
public:
    /**
     * @param num Number of objective vectors.
     * @param[out] store Pointer to an output array of length `num`.
     * After `finish()` is called, this will contain the minimum/maximum for each objective vector.
     * @param skip_nan See `Options::skip_nan` for details.
     * @param subtract Offset to subtract from each element of `index` before using it to index into `store`.
     * Only relevant if `store` holds statistics for a contiguous subset of objective vectors,
     * e.g., during task allocation for parallelization.
     */
    RunningSparse(Index_ num, Output_* store, bool skip_nan, Index_ subtract = 0) : 
        my_num(num), my_store(store), my_skip_nan(skip_nan), my_subtract(subtract) {}

    /**
     * Add the next observed vector to the min/max calculation.
     * @param[in] value Value of structural non-zero elements.
     * @param[in] index Index of structural non-zero elements.
     * @param number Number of non-zero elements in `value` and `index`.
     */
    void add(const Value_* value, const Index_* index, Index_ number) {
        if (my_count == 0) {
            my_nonzero.resize(my_num);
            std::fill_n(my_store, my_num, internal::choose_placeholder<minimum_, Value_>());

            if (!my_skip_nan) {
                for (Index_ i = 0; i < number; ++i, ++value, ++index) {
                    auto val = *value;
                    auto idx = *index - my_subtract;
                    my_store[idx] = val;
                    ++my_nonzero[idx];
                }
                my_count = 1;
                return;
            }
        }

        for (Index_ i = 0; i < number; ++i, ++value, ++index) {
            auto val = *value;
            auto idx = *index - my_subtract;
            auto& current = my_store[idx];
            if (internal::is_better<minimum_>(current, val)) { // this should implicitly skip NaNs, any NaN comparison will be false.
                current = val;
            }
            ++my_nonzero[idx];
        }

        ++my_count;
    }

    /**
     * Finish the min/max calculation once all observed vectors have been passed to `add()`. 
     */
    void finish() {
        if (my_count) {
            for (Index_ i = 0; i < my_num; ++i) {
                if (my_count > my_nonzero[i]) {
                    auto& current = my_store[i];
                    if (internal::is_better<minimum_>(current, 0)) {
                        current = 0;
                    }
                }
            }
        } else {
            std::fill_n(my_store, my_num, internal::choose_placeholder<minimum_, Value_>());
        }
    }

private:
    Index_ my_num;
    Output_* my_store;
    bool my_skip_nan;
    Index_ my_subtract;
    Index_ my_count = 0;
    std::vector<Index_> my_nonzero;
};

/**
 * Compute ranges for each element of a chosen dimension of a `tatami::Matrix`.
 *
 * @tparam Value_ Type of the matrix value, should be numeric.
 * @tparam Index_ Type of the row/column indices.
 * @tparam Output_ Type of the output value.
 *
 * @param row Whether to compute variances for the rows.
 * @param p Pointer to a `tatami::Matrix`.
 * @param[out] min_out Pointer to an array of length equal to the number of rows (if `row = true`) or columns (otherwise).
 * On output, this will contain the minimum of each row/column.
 * Alternatively, this may be NULL, in which case the minima are not computed.
 * @param[out] max_out Pointer to an array of length equal to the number of rows (if `row = true`) or columns (otherwise).
 * On output, this will contain the maximum of each row/column.
 * Alternatively, this may be NULL, in which case the maxima are not computed.
 * @param ropt Range calculation options.
 */
template<typename Value_, typename Index_, typename Output_>
void apply(bool row, const tatami::Matrix<Value_, Index_>* p, Output_* min_out, Output_* max_out, const Options& ropt) {
    auto dim = (row ? p->nrow() : p->ncol());
    auto otherdim = (row ? p->ncol() : p->nrow());
    const bool direct = p->prefer_rows() == row;

    bool store_min = min_out != NULL;
    bool store_max = max_out != NULL;

    if (p->sparse()) {
        tatami::Options opt;
        opt.sparse_ordered_index = false;

        if (direct) {
            opt.sparse_extract_index = false;
            tatami::parallelize([&](size_t, Index_ s, Index_ l) {
                auto ext = tatami::consecutive_extractor<true>(p, row, s, l, opt);
                std::vector<Value_> vbuffer(otherdim);
                for (Index_ x = 0; x < l; ++x) {
                    auto out = ext->fetch(vbuffer.data(), NULL);
                    if (store_min) {
                        min_out[x + s] = ranges::direct<true>(out.value, out.number, otherdim, ropt.skip_nan);
                    }
                    if (store_max) {
                        max_out[x + s] = ranges::direct<false>(out.value, out.number, otherdim, ropt.skip_nan);
                    }
                }
            }, dim, ropt.num_threads);

        } else {
            tatami::parallelize([&](size_t thread, Index_ s, Index_ l) {
                auto ext = tatami::consecutive_extractor<true>(p, !row, 0, otherdim, s, l, opt);
                std::vector<Value_> vbuffer(l);
                std::vector<Index_> ibuffer(l);

                auto local_min = (store_min ? LocalOutputBuffer<Output_>(thread, s, l, min_out) : LocalOutputBuffer<Output_>());
                auto local_max = (store_max ? LocalOutputBuffer<Output_>(thread, s, l, max_out) : LocalOutputBuffer<Output_>());
                ranges::RunningSparse<true, Output_, Value_, Index_> runmin(l, local_min.data(), ropt.skip_nan, s);
                ranges::RunningSparse<false, Output_, Value_, Index_> runmax(l, local_max.data(), ropt.skip_nan, s);

                for (Index_ x = 0; x < otherdim; ++x) {
                    auto out = ext->fetch(vbuffer.data(), ibuffer.data());
                    if (store_min) {
                        runmin.add(out.value, out.index, out.number);
                    }
                    if (store_max) {
                        runmax.add(out.value, out.index, out.number);
                    }
                }

                if (store_min) {
                    runmin.finish();
                    local_min.transfer();
                }
                if (store_max) {
                    runmax.finish();
                    local_max.transfer();
                }
            }, dim, ropt.num_threads);
        }

    } else {
        if (direct) {
            tatami::parallelize([&](size_t, Index_ s, Index_ l) {
                auto ext = tatami::consecutive_extractor<false>(p, row, s, l);
                std::vector<Value_> buffer(otherdim);
                for (Index_ x = 0; x < l; ++x) {
                    auto ptr = ext->fetch(buffer.data());
                    if (store_min) {
                        min_out[x + s] = ranges::direct<true>(ptr, otherdim, ropt.skip_nan);
                    }
                    if (store_max) {
                        max_out[x + s] = ranges::direct<false>(ptr, otherdim, ropt.skip_nan);
                    }
                }
            }, dim, ropt.num_threads);

        } else {
            tatami::parallelize([&](size_t thread, Index_ s, Index_ l) {
                auto ext = tatami::consecutive_extractor<false>(p, !row, 0, otherdim, s, l);
                std::vector<Value_> buffer(l);

                auto local_min = (store_min ? LocalOutputBuffer<Output_>(thread, s, l, min_out) : LocalOutputBuffer<Output_>());
                auto local_max = (store_max ? LocalOutputBuffer<Output_>(thread, s, l, max_out) : LocalOutputBuffer<Output_>());
                ranges::RunningDense<true, Output_, Value_, Index_> runmin(l, local_min.data(), ropt.skip_nan);
                ranges::RunningDense<false, Output_, Value_, Index_> runmax(l, local_max.data(), ropt.skip_nan);

                for (Index_ x = 0; x < otherdim; ++x) {
                    auto ptr = ext->fetch(buffer.data());
                    if (store_min) {
                        runmin.add(ptr);
                    }
                    if (store_max) {
                        runmax.add(ptr);
                    }
                }

                if (store_min) {
                    runmin.finish();
                    local_min.transfer();
                }
                if (store_max) {
                    runmax.finish();
                    local_max.transfer();
                }
            }, dim, ropt.num_threads);
        }
    }

    return;
}

/**
 * Wrapper around `apply()` for column ranges.
 *
 * @tparam Output Type of the output value.
 * @tparam Value_ Type of the matrix value.
 * @tparam Index_ Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param ropt Range calculation options.
 *
 * @return A pair of vectors, each of length equal to the number of columns.
 * The first and second vector contains the minimum and maximum value per column, respectively.
 */
template<typename Output_ = double, typename Value_, typename Index_>
std::pair<std::vector<Output_>, std::vector<Output_> > by_column(const tatami::Matrix<Value_, Index_>* p, const Options& ropt) {
    std::vector<Output_> mins(p->ncol()), maxs(p->ncol());
    apply(false, p, mins.data(), maxs.data(), ropt);
    return std::make_pair(std::move(mins), std::move(maxs));
}

/**
 * Overload with default options.
 *
 * @tparam Output Type of the output value.
 * @tparam Value_ Type of the matrix value.
 * @tparam Index_ Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::Matrix`.
 *
 * @return A pair of vectors, each of length equal to the number of columns.
 * The first and second vector contains the minimum and maximum value per column, respectively.
 */
template<typename Output_ = double, typename Value_, typename Index_>
std::pair<std::vector<Output_>, std::vector<Output_> > by_column(const tatami::Matrix<Value_, Index_>* p) {
    return by_column<Output_>(p, Options());
}

/**
 * Wrapper around `apply()` for row ranges.
 *
 * @tparam Output Type of the output value.
 * @tparam Value_ Type of the matrix value.
 * @tparam Index_ Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param ropt Range calculation options.
 *
 * @return A pair of vectors, each of length equal to the number of rows.
 * The first and second vector contains the minimum and maximum value per row, respectively.
 */
template<typename Output_ = double, typename Value_, typename Index_>
std::pair<std::vector<Output_>, std::vector<Output_> > by_row(const tatami::Matrix<Value_, Index_>* p, const Options& ropt) {
    std::vector<Output_> mins(p->nrow()), maxs(p->nrow());
    apply(true, p, mins.data(), maxs.data(), ropt);
    return std::make_pair(std::move(mins), std::move(maxs));
}

/**
 * Overload with default options.
 *
 * @tparam Output Type of the output value.
 * @tparam Value_ Type of the matrix value.
 * @tparam Index_ Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::Matrix`.
 *
 * @return A pair of vectors, each of length equal to the number of rows.
 * The first and second vector contains the minimum and maximum value per row, respectively.
 */
template<typename Output_ = double, typename Value_, typename Index_>
std::pair<std::vector<Output_>, std::vector<Output_> > by_row(const tatami::Matrix<Value_, Index_>* p) {
    return by_row<Output_>(p, Options());
}

}

}

#endif
