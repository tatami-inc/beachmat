#ifndef TATAMI_STATS__SUMS_HPP
#define TATAMI_STATS__SUMS_HPP

#include "tatami/tatami.hpp"
#include "utils.hpp"

#include <vector>
#include <numeric>
#include <algorithm>

/**
 * @file sums.hpp
 *
 * @brief Compute row and column sums from a `tatami::Matrix`.
 */

namespace tatami_stats {

/**
 * @brief Functions for computing dimension-wise sums.
 * @namespace tatami_stats::sums
 */
namespace sums {

/**
 * @brief Summation options.
 */
struct Options {
    /**
     * Whether to check for NaNs in the input, and skip them.
     * If false, NaNs are assumed to be absent, and the behavior of summation in the presence of NaNs is undefined.
     */
    bool skip_nan = false;

    /**
     * Number of threads to use when computing sums across a `tatami::Matrix`.
     */
    int num_threads = 1;
};

/**
 * Directly sum an array of values using naive accumulation.
 * This is best used with a sufficiently high-precision `Output_`, hence the default of `double`.
 *
 * @tparam Output_ Type of the output data.
 * @tparam Value_ Type of the input data.
 * @tparam Index_ Type of the row/column index.
 *
 * @param[in] ptr Pointer to an array of values of length `num`.
 * @param num Size of the array.
 * @param skip_nan See `Options::skip_nan`.
 * @return The sum.
 */
template<typename Output_ = double, typename Value_, typename Index_>
Output_ direct(const Value_* ptr, Index_ num, bool skip_nan) {
    if (skip_nan) {
        Output_ sum = 0;
        for (Index_ i = 0; i < num; ++i) {
            auto val = ptr[i];
            if (!std::isnan(val)) {
                sum += val;
            }
        }
        return sum;
    } else {
        return std::accumulate(ptr, ptr + num, static_cast<Output_>(0));
    }
}

/**
 * @brief Running sums from dense data.
 *
 * This considers a scenario with a set of equilength "objective" vectors [V1, V2, V3, ..., Vn],
 * but data are only available for "observed" vectors [P1, P2, P3, ..., Pm],
 * where Pi[j] contains the i-th element of objective vector Vj.
 * The idea is to repeatedly call `add()` for `ptr` corresponding to observed vectors from 0 to m - 1,
 * and then finally call `finish()` to obtain the sum for each objective vector.
 *
 * This class uses naive accumulation to obtain the sum for each objective vector.
 * Callers should use a sufficiently high-precision `Output_` such as `double` to mitigate round-off errors.
 *
 * @tparam Output_ Type of the output data.
 * @tparam Value_ Type of the input value.
 * @tparam Index_ Type of the row/column indices.
 */
template<typename Output_, typename Value_, typename Index_>
class RunningDense {
public:
    /**
     * @param num Number of objective vectors, i.e., n.
     * @param[out] sum Pointer to an output array of length `num`.
     * This should be zeroed on input, and will store the running sums after each `add()`.
     * @param skip_nan See `Options::skip_nan` for details.
     */
    RunningDense(Index_ num, Output_* sum, bool skip_nan) : my_num(num), my_sum(sum), my_skip_nan(skip_nan) {}

    /**
     * Add the next observed vector to the running sums.
     * @param[in] ptr Pointer to an array of values of length `my_num`, corresponding to an observed vector.
     */
    void add(const Value_* ptr) {
        if (my_skip_nan) {
            for (Index_ i = 0; i < my_num; ++i) {
                auto val = ptr[i];
                if (!std::isnan(val)) {
                    my_sum[i] += val;
                }
            }
        } else {
            for (Index_ i = 0; i < my_num; ++i) {
                my_sum[i] += ptr[i];
            }
        }
    }

private:
    Index_ my_num;
    Output_* my_sum;
    bool my_skip_nan;
};

/**
 * @brief Running sums from sparse data.
 *
 * Compute running sums from sparse data. 
 * This is the counterpart to `RunningDense`, but for sparse observed vectors.
 *
 * @tparam Output_ Type of the output data.
 * @tparam Value_ Type of the input value.
 * @tparam Index_ Type of the row/column indices.
 */
template<typename Output_, typename Value_, typename Index_>
class RunningSparse {
public:
    /**
     * @param[out] sum Pointer to an output array of length equal to the number of objective vectors.
     * This should be zeroed on input, and will store the running sums after each `add()`.
     * @param skip_nan See `Options::skip_nan` for details.
     * @param subtract Offset to subtract from each element of `index` before using it to index into `mean` and friends.
     * Only relevant if `mean` and friends hold statistics for a contiguous subset of objective vectors,
     * e.g., during task allocation for parallelization.
     */
    RunningSparse(Output_* sum, bool skip_nan, Index_ subtract = 0) : 
        my_sum(sum), my_skip_nan(skip_nan), my_subtract(subtract) {}

    /**
     * Add the next observed vector to the running sums.
     * @param[in] value Value of structural non-zero elements.
     * @param[in] index Index of structural non-zero elements.
     * This does not have to be sorted.
     * @param number Number of non-zero elements in `value` and `index`.
     */
    void add(const Value_* value, const Index_* index, Index_ number) {
        if (my_skip_nan) {
            for (Index_ i = 0; i < number; ++i) {
                auto val = value[i];
                if (!std::isnan(val)) {
                    my_sum[index[i] - my_subtract] += val;
                }
            }
        } else {
            for (Index_ i = 0; i < number; ++i) {
                my_sum[index[i] - my_subtract] += value[i];
            }
        }
    }

private:
    Output_* my_sum;
    bool my_skip_nan;
    Index_ my_subtract;
};

/**
 * Compute sums for each element of a chosen dimension of a `tatami::Matrix`.
 * This is done using `direct()`, `RunningDense` or `RunningSparse`, 
 * depending on the requested dimension in `row` and the preferred dimension for access in `p`.
 * Note that all sums are obtained using naive accumulation,
 * so it is best to use a sufficiently high-precision `Output_` to mitigate round-off errors.
 *
 * @tparam Value_ Type of the matrix value, should be numeric.
 * @tparam Index_ Type of the row/column indices.
 * @tparam Output_ Type of the output value.
 *
 * @param row Whether to compute variances for the rows.
 * @param p Pointer to a `tatami::Matrix`.
 * @param[out] output Pointer to an array of length equal to the number of rows (if `row = true`) or columns (otherwise).
 * On output, this will contain the row/column variances.
 * @param sopt Summation options.
 */
template<typename Value_, typename Index_, typename Output_>
void apply(bool row, const tatami::Matrix<Value_, Index_>* p, Output_* output, const Options& sopt) {
    auto dim = (row ? p->nrow() : p->ncol());
    auto otherdim = (row ? p->ncol() : p->nrow());
    const bool direct = p->prefer_rows() == row;

    if (p->sparse()) {
        if (direct) {
            tatami::Options opt;
            opt.sparse_extract_index = false;

            tatami::parallelize([&](size_t, Index_ s, Index_ l) {
                auto ext = tatami::consecutive_extractor<true>(p, row, s, l, opt);
                std::vector<Value_> vbuffer(otherdim);
                for (Index_ x = 0; x < l; ++x) {
                    auto out = ext->fetch(vbuffer.data(), NULL);
                    output[x + s] = sums::direct(out.value, out.number, sopt.skip_nan);
                }
            }, dim, sopt.num_threads);

        } else {
            tatami::Options opt;
            opt.sparse_ordered_index = false;

            tatami::parallelize([&](size_t thread, Index_ s, Index_ l) {
                auto ext = tatami::consecutive_extractor<true>(p, !row, 0, otherdim, s, l, opt);
                std::vector<Value_> vbuffer(l);
                std::vector<Index_> ibuffer(l);

                LocalOutputBuffer<Output_> local_output(thread, s, l, output);
                sums::RunningSparse<Output_, Value_, Index_> runner(local_output.data(), sopt.skip_nan, s);

                for (Index_ x = 0; x < otherdim; ++x) {
                    auto out = ext->fetch(vbuffer.data(), ibuffer.data());
                    runner.add(out.value, out.index, out.number);
                }

                local_output.transfer();
            }, dim, sopt.num_threads);
        }

    } else {
        if (direct) {
            tatami::parallelize([&](size_t, Index_ s, Index_ l) {
                auto ext = tatami::consecutive_extractor<false>(p, row, s, l);
                std::vector<Value_> buffer(otherdim);
                for (Index_ x = 0; x < l; ++x) {
                    auto out = ext->fetch(buffer.data());
                    output[x + s] = sums::direct(out, otherdim, sopt.skip_nan);
                }
            }, dim, sopt.num_threads);

        } else {
            tatami::parallelize([&](size_t thread, Index_ s, Index_ l) {
                auto ext = tatami::consecutive_extractor<false>(p, !row, 0, otherdim, s, l);
                std::vector<Value_> buffer(l);

                LocalOutputBuffer<Output_> local_output(thread, s, l, output);
                sums::RunningDense<Output_, Value_, Index_> runner(l, local_output.data(), sopt.skip_nan);

                for (Index_ x = 0; x < otherdim; ++x) {
                    auto out = ext->fetch(buffer.data());
                    runner.add(out);
                }

                local_output.transfer();
            }, dim, sopt.num_threads);
        }
    }

    return;
}

/**
 * Wrapper around `apply()` for column sums.
 *
 * @tparam Value_ Type of the matrix value, should be summable.
 * @tparam Index_ Type of the row/column indices.
 * @tparam Output_ Type of the output value.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param sopt Summation options.
 * @return Vector containing the sum for each column.
 */
template<typename Output_ = double, typename Value_, typename Index_>
std::vector<Output_> by_column(const tatami::Matrix<Value_, Index_>* p, const Options& sopt) {
    std::vector<Output_> output(p->ncol());
    apply(false, p, output.data(), sopt);
    return output;
}

/**
 * Overload with default options.
 *
 * @tparam Value_ Type of the matrix value, should be summable.
 * @tparam Index_ Type of the row/column indices.
 * @tparam Output_ Type of the output value.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @return Vector containing the sum for each column.
 */
template<typename Output_ = double, typename Value_, typename Index_>
std::vector<Output_> by_column(const tatami::Matrix<Value_, Index_>* p) {
    return by_column(p, Options());
}

/**
 * Wrapper around `apply()` for row sums.
 *
 * @tparam Output_ Type of the output value.
 * @tparam Value_ Type of the matrix value, should be summable.
 * @tparam Index_ Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param sopt Summation options.
 * @return Vector containing the sum of each row.
 */
template<typename Output_ = double, typename Value_, typename Index_>
std::vector<Output_> by_row(const tatami::Matrix<Value_, Index_>* p, const Options& sopt) {
    std::vector<Output_> output(p->nrow());
    apply(true, p, output.data(), sopt);
    return output;
}

/**
 * Overload with default options.
 *
 * @tparam Output_ Type of the output value.
 * @tparam Value_ Type of the matrix value, should be summable.
 * @tparam Index_ Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @return Vector containing the sum of each row.
 */
template<typename Output_ = double, typename Value_, typename Index_>
std::vector<Output_> by_row(const tatami::Matrix<Value_, Index_>* p) {
    return by_row(p, Options());
}

}

}

#endif
