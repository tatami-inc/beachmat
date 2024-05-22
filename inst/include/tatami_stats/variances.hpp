#ifndef TATAMI_STATS_VARS_HPP
#define TATAMI_STATS_VARS_HPP

#include "tatami/tatami.hpp"
#include "utils.hpp"

#include <vector>
#include <cmath>
#include <numeric>
#include <limits>
#include <algorithm>

/**
 * @file variances.hpp
 *
 * @brief Compute row and column variances from a `tatami::Matrix`.
 */

namespace tatami_stats {

/**
 * @brief Functions for computing dimension-wise variances.
 * @namespace tatami_stats::variances
 */
namespace variances {

/**
 * @brief Variance calculation options.
 */
struct Options {
    /**
     * Whether to check for NaNs in the input, and skip them.
     * If false, NaNs are assumed to be absent, and the behavior of the variance calculation in the presence of NaNs is undefined.
     */
    bool skip_nan = false;

    /**
     * Number of threads to use when computing variances across a `tatami::Matrix`.
     */
    int num_threads = 1;
};

/**
 * @cond
 */
namespace internal {

template<typename Output_ = double, typename Value_, typename Index_ >
void add_welford(Output_& mean, Output_& sumsq, Value_ value, Index_ count) {
    Output_ delta = value - mean;
    mean += delta / count;
    sumsq += delta * (value - mean);
}

template<typename Output_ = double, typename Index_ >
void add_welford_zeros(Output_& mean, Output_& sumsq, Index_ num_nonzero, Index_ num_all) {
    auto ratio = static_cast<Output_>(num_nonzero) / static_cast<Output_>(num_all);
    sumsq += mean * mean * ratio * (num_all - num_nonzero);
    mean *= ratio;
}

}
/**
 * @endcond
 */

/**
 * Compute the mean and variance from a sparse array of values.
 * This uses the standard two-pass algorithm with naive accumulation of the sum of squared differences;
 * thus, it is best used with a sufficiently high-precision `Output_` like `double`.
 *
 * @tparam Output_ Type of the output data.
 * @tparam Value_ Type of the input data.
 * @tparam Index_ Type of the row/column indices.
 *
 * @param[in] value Pointer to an array of values of length `num`.
 * @param num_nonzero Length of the array pointed to by `value`.
 * @param num_all Total number of values in the dataset, including the zeros not in `value`.
 * This should be greater than or equal to `num_nonzero`.
 * @param skip_nan See `Options::skip_nan`.
 *
 * @return The sample mean and variance of values in the sparse array.
 * This may be NaN if there are not enough (non-NaN) values in `value`.
 */
template<typename Output_ = double, typename Value_, typename Index_ >
std::pair<Output_, Output_> direct(const Value_* value, Index_ num_nonzero, Index_ num_all, bool skip_nan) {
    Output_ mean = 0;
    Index_ lost = 0;

    if (skip_nan) {
        auto copy = value;
        for (Index_ i = 0; i < num_nonzero; ++i, ++copy) {
            auto val = *copy;
            if (std::isnan(val)) {
                ++lost;
            } else {
                mean += val;
            }
        }
    } else {
        auto copy = value;
        for (Index_ i = 0; i < num_nonzero; ++i, ++copy) {
            mean += *copy;
        }
    }

    auto count = num_all - lost;
    mean /= count;

    Output_ var = 0;
    if (skip_nan) {
        for (Index_ i = 0; i < num_nonzero; ++i, ++value) {
            auto val = *value;
            if (!std::isnan(val)) {
                var += (val - mean) * (val - mean);
            }
        }
    } else {
        for (Index_ i = 0; i < num_nonzero; ++i, ++value) {
            auto val = *value;
            var += (val - mean) * (val - mean);
        }
    }

    if (num_nonzero < num_all) {
        var += (num_all - num_nonzero) * mean * mean;
    }

    if (count == 0) {
        return std::make_pair(std::numeric_limits<Output_>::quiet_NaN(), std::numeric_limits<Output_>::quiet_NaN());
    } else if (count == 1) {
        return std::make_pair(mean, std::numeric_limits<Output_>::quiet_NaN());
    } else {
        return std::make_pair(mean, var / (count - 1));
    }
}

/**
 * Compute the mean and variance from an array of values.
 * This uses the standard two-pass algorithm with naive accumulation of the sum of squared differences;
 * thus, it is best used with a sufficiently high-precision `Output_` like `double`.
 *
 * @tparam Output_ Type of the output data.
 * @tparam Value_ Type of the input data.
 * @tparam Index_ Type of the row/column indices.
 *
 * @param[in] ptr Pointer to an array of values of length `num`.
 * @param num Size of the array.
 * @param skip_nan See `Options::skip_nan`.
 *
 * @return The sample mean and variance of values in `[ptr, ptr + num)`.
 * This may be NaN if there are not enough (non-NaN) values in `ptr`.
 */
template<typename Output_ = double, typename Value_, typename Index_ >
std::pair<Output_, Output_> direct(const Value_* ptr, Index_ num, bool skip_nan) {
    return direct<Output_>(ptr, num, num, skip_nan);
}

/**
 * @brief Running variances from dense data.
 *
 * Compute running means and variances from dense data using Welford's method.
 * This considers a scenario with a set of equilength "objective" vectors [V1, V2, V3, ..., Vn],
 * but data are only available for "observed" vectors [P1, P2, P3, ..., Pm],
 * where Pi[j] contains the i-th element of objective vector Vj.
 * The idea is to repeatedly call `add()` for `ptr` corresponding to observed vectors from 0 to m - 1,
 * and then finally call `finish()` to obtain the mean and variance for each objective vector.
 *
 * @tparam Output_ Type of the output data.
 * @tparam Value_ Type of the input data.
 * @tparam Index_ Type of the row/column indices.
 */
template<typename Output_, typename Value_, typename Index_>
class RunningDense {
public:
    /**
     * @param num Number of objective vectors, i.e., n.
     * @param[out] mean Pointer to an output array of length `num`.
     * This should be zeroed on input; after `finish()` is called, this will contain the means for each objective vector.
     * @param[out] variance Pointer to an output array of length `num`, containing the variances for each objective vector.
     * This should be zeroed on input; after `finish()` is called, this will contain the sample variance for each objective vector.
     * @param skip_nan See `Options::skip_nan` for details.
     */
    RunningDense(Index_ num, Output_* mean, Output_* variance, bool skip_nan) : 
        my_num(num), my_mean(mean), my_variance(variance), my_skip_nan(skip_nan), my_ok_count(skip_nan ? num : 0) {}

    /**
     * Add the next observed vector to the variance calculation.
     * @param[in] ptr Pointer to an array of values of length `num`, corresponding to an observed vector.
     */
    void add(const Value_* ptr) {
        if (my_skip_nan) {
            for (Index_ i = 0; i < my_num; ++i, ++ptr) {
                auto val = *ptr;
                if (!std::isnan(val)) {
                    internal::add_welford(my_mean[i], my_variance[i], val, ++(my_ok_count[i]));
                }
            }
        } else {
            ++my_count;
            for (Index_ i = 0; i < my_num; ++i, ++ptr) {
                internal::add_welford(my_mean[i], my_variance[i], *ptr, my_count);
            }
        }
    }

    /**
     * Finish the variance calculation once all observed vectors have been passed to `add()`. 
     */
    void finish() {
        if (my_skip_nan) {
            for (Index_ i = 0; i < my_num; ++i) {
                auto ct = my_ok_count[i];
                if (ct < 2) {
                    my_variance[i] = std::numeric_limits<Output_>::quiet_NaN();
                    if (ct == 0) {
                        my_mean[i] = std::numeric_limits<Output_>::quiet_NaN();
                    }
                } else {
                    my_variance[i] /= ct - 1;
                }
            }
        } else {
            if (my_count < 2) {
                std::fill_n(my_variance, my_num, std::numeric_limits<Output_>::quiet_NaN());
                if (my_count == 0) {
                    std::fill_n(my_mean, my_num, std::numeric_limits<Output_>::quiet_NaN());
                }
            } else {
                for (Index_ i = 0; i < my_num; ++i) {
                    my_variance[i] /= my_count - 1;
                }
            }
        }
    }

private:
    Index_ my_num;
    Output_* my_mean;
    Output_* my_variance;
    bool my_skip_nan;
    Index_ my_count = 0;
    std::vector<Index_> my_ok_count;
};

/**
 * @brief Running variances from sparse data.
 *
 * Compute running means and variances from sparse data using Welford's method.
 * This does the same as its dense overload for sparse observed vectors.
 *
 * @tparam Output_ Type of the output data.
 * @tparam Value_ Type of the input data.
 * @tparam Index_ Type of the row/column indices.
 */
template<typename Output_, typename Value_, typename Index_>
class RunningSparse {
public:
    /**
     * @param num Number of objective vectors.
     * @param[out] mean Pointer to an output array of length `num`, containing the means for each objective vector.
     * This should be zeroed on input; after `finish()` is called, this will contain the mean for each objective vector.
     * @param[out] variance Pointer to an output array of length `num`, containing the variances for each objective vector.
     * This should be zeroed on input; after `finish()` is called, this will contain the sample variance for each objective vector.
     * @param skip_nan See `Options::skip_nan` for details.
     * @param subtract Offset to subtract from each element of `index` before using it to index into `mean` and friends.
     * Only relevant if `mean` and friends hold statistics for a contiguous subset of objective vectors,
     * e.g., during task allocation for parallelization.
     */
    RunningSparse(Index_ num, Output_* mean, Output_* variance, bool skip_nan, Index_ subtract = 0) : 
        my_num(num), my_mean(mean), my_variance(variance), my_nonzero(num), my_skip_nan(skip_nan), my_subtract(subtract), my_nan(skip_nan ? num : 0) {}

    /**
     * Add the next observed vector to the variance calculation.
     * @param[in] value Value of structural non-zero elements.
     * @param[in] index Index of structural non-zero elements.
     * @param number Number of non-zero elements in `value` and `index`.
     */
    void add(const Value_* value, const Index_* index, Index_ number) {
        ++my_count;
        if (my_skip_nan) {
            for (Index_ i = 0; i < number; ++i) {
                auto val = value[i];
                auto ri = index[i] - my_subtract;
                if (std::isnan(val)) {
                    ++my_nan[ri];
                } else {
                    internal::add_welford(my_mean[ri], my_variance[ri], val, ++(my_nonzero[ri]));
                }
            }

        } else {
            for (Index_ i = 0; i < number; ++i) {
                auto ri = index[i] - my_subtract;
                internal::add_welford(my_mean[ri], my_variance[ri], value[i], ++(my_nonzero[ri]));
            }
        }
    }

    /**
     * Finish the variance calculation once all observed vectors have been passed to `add()`. 
     */
    void finish() {
        if (my_skip_nan) {
            for (Index_ i = 0; i < my_num; ++i) {
                auto& curM = my_mean[i];
                auto& curV = my_variance[i];
                auto ct = my_count - my_nan[i];

                if (ct < 2) {
                    curV = std::numeric_limits<Output_>::quiet_NaN();
                    if (ct == 0) {
                        curM = std::numeric_limits<Output_>::quiet_NaN();
                    }
                } else {
                    internal::add_welford_zeros(curM, curV, my_nonzero[i], ct);
                    curV /= ct - 1;
                }
            }

        } else {
            if (my_count < 2) {
                std::fill_n(my_variance, my_num, std::numeric_limits<Output_>::quiet_NaN());
                if (my_count == 0) {
                    std::fill_n(my_mean, my_num, std::numeric_limits<Output_>::quiet_NaN());
                }
            } else {
                for (Index_ i = 0; i < my_num; ++i) {
                    auto& var = my_variance[i];
                    internal::add_welford_zeros(my_mean[i], var, my_nonzero[i], my_count);
                    var /= my_count - 1;
                }
            }
        }
    }

private:
    Index_ my_num;
    Output_* my_mean;
    Output_* my_variance;
    std::vector<Index_> my_nonzero;
    bool my_skip_nan;
    Index_ my_subtract;
    Index_ my_count = 0;
    std::vector<Index_> my_nan;
};

/**
 * Compute variances for each element of a chosen dimension of a `tatami::Matrix`.
 * This may use either Welford's method or the standard two-pass method,
 * depending on the dimension in `row` and the preferred access dimension of `p`.
 *
 * @tparam Value_ Type of the matrix value, should be numeric.
 * @tparam Index_ Type of the row/column indices.
 * @tparam Output_ Type of the output value.
 *
 * @param row Whether to compute variances for the rows.
 * @param p Pointer to a `tatami::Matrix`.
 * @param[out] output Pointer to an array of length equal to the number of rows (if `row = true`) or columns (otherwise).
 * On output, this will contain the row/column variances.
 * @param vopt Variance calculation options.
 */
template<typename Value_, typename Index_, typename Output_>
void apply(bool row, const tatami::Matrix<Value_, Index_>* p, Output_* output, const Options& vopt) {
    auto dim = (row ? p->nrow() : p->ncol());
    auto otherdim = (row ? p->ncol() : p->nrow());
    const bool direct = p->prefer_rows() == row;

    if (p->sparse()) {
        if (direct) {
            tatami::Options opt;
            opt.sparse_extract_index = false;
            tatami::parallelize([&](size_t, Index_ s, Index_ l) {
                auto ext = tatami::consecutive_extractor<true>(p, row, s, l);
                std::vector<Value_> vbuffer(otherdim);
                for (Index_ x = 0; x < l; ++x) {
                    auto out = ext->fetch(vbuffer.data(), NULL);
                    output[x + s] = variances::direct<Output_>(out.value, out.number, otherdim, vopt.skip_nan).second;
                }
            }, dim, vopt.num_threads);

        } else {
            tatami::parallelize([&](size_t thread, Index_ s, Index_ l) {
                auto ext = tatami::consecutive_extractor<true>(p, !row, 0, otherdim, s, l);
                std::vector<Value_> vbuffer(l);
                std::vector<Index_> ibuffer(l);

                std::vector<Output_> running_means(l);
                LocalOutputBuffer<Output_> local_output(thread, s, l, output);
                variances::RunningSparse<Output_, Value_, Index_> runner(l, running_means.data(), local_output.data(), vopt.skip_nan, s);

                for (Index_ x = 0; x < otherdim; ++x) {
                    auto out = ext->fetch(vbuffer.data(), ibuffer.data());
                    runner.add(out.value, out.index, out.number);
                }
                runner.finish();

                local_output.transfer(); 
            }, dim, vopt.num_threads);
        }

    } else {
        if (direct) {
            tatami::parallelize([&](size_t, Index_ s, Index_ l) {
                auto ext = tatami::consecutive_extractor<false>(p, row, s, l);
                std::vector<Value_> buffer(otherdim);
                for (Index_ x = 0; x < l; ++x) {
                    auto out = ext->fetch(buffer.data());
                    output[x + s] = variances::direct<Output_>(out, otherdim, vopt.skip_nan).second;
                }
            }, dim, vopt.num_threads);

        } else {
            tatami::parallelize([&](size_t thread, Index_ s, Index_ l) {
                auto ext = tatami::consecutive_extractor<false>(p, !row, 0, otherdim, s, l);
                std::vector<Value_> buffer(l);

                std::vector<Output_> running_means(l);
                LocalOutputBuffer<Output_> local_output(thread, s, l, output);
                variances::RunningDense<Output_, Value_, Index_> runner(l, running_means.data(), local_output.data(), vopt.skip_nan);

                for (Index_ x = 0; x < otherdim; ++x) {
                    runner.add(ext->fetch(buffer.data()));
                }
                runner.finish();

                local_output.transfer(); 
            }, dim, vopt.num_threads);
        }
    }
}

/**
 * Wrapper around `apply()` for column variances.
 *
 * @tparam Output_ Type of the output value.
 * @tparam Value_ Type of the matrix value, should be numeric.
 * @tparam Index_ Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param vopt Variance calculation options.
 *
 * @return A vector of length equal to the number of columns, containing the column variances.
 */
template<typename Output_ = double, typename Value_, typename Index_>
std::vector<Output_> by_column(const tatami::Matrix<Value_, Index_>* p, const Options& vopt) {
    std::vector<Output_> output(p->ncol());
    apply(false, p, output.data(), vopt);
    return output;
}

/**
 * Overload with default options.
 *
 * @tparam Output_ Type of the output value.
 * @tparam Value_ Type of the matrix value, should be numeric.
 * @tparam Index_ Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::Matrix`.
 *
 * @return A vector of length equal to the number of columns, containing the column variances.
 */
template<typename Output_ = double, typename Value_, typename Index_>
std::vector<Output_> by_column(const tatami::Matrix<Value_, Index_>* p) {
    Options vopt;
    return by_column(p, vopt);
}

/**
 * Wrapper around `apply()` for column variances.
 *
 * @tparam Output_ Type of the output value.
 * @tparam Value_ Type of the matrix value, should be numeric.
 * @tparam Index_ Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::Matrix`.
 * @param vopt Variance calculation options.
 *
 * @return A vector of length equal to the number of rows, containing the row variances.
 */
template<typename Output_ = double, typename Value_, typename Index_>
std::vector<Output_> by_row(const tatami::Matrix<Value_, Index_>* p, const Options& vopt) {
    std::vector<Output_> output(p->nrow());
    apply(true, p, output.data(), vopt);
    return output;
}

/**
 * Overload with default options.
 *
 * @tparam Output_ Type of the output value.
 * @tparam Value_ Type of the matrix value, should be numeric.
 * @tparam Index_ Type of the row/column indices.
 *
 * @param p Pointer to a `tatami::Matrix`.
 *
 * @return A vector of length equal to the number of rows, containing the row variances.
 */
template<typename Output_ = double, typename Value_, typename Index_>
std::vector<Output_> by_row(const tatami::Matrix<Value_, Index_>* p) {
    Options vopt;
    return by_row(p, vopt);
}

}

}

#endif
