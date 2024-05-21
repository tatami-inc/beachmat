#ifndef TATAMI_PARALLELIZE_HPP
#define TATAMI_PARALLELIZE_HPP

#include <vector>
#include <cmath>

#ifndef TATAMI_CUSTOM_PARALLEL
#ifndef _OPENMP
#include <thread>
#endif
#include <string>
#include <stdexcept>
#endif

/**
 * @file parallelize.hpp
 *
 * @brief Parallelized iteration over a `tatami::Matrix`.
 */

namespace tatami {

/**
 * Apply a function to a set of tasks in parallel, usually for iterating over a `Matrix`.
 * This can be done using:
 *
 * - OpenMP, if available and enabled by the compiler.
 * - Using a custom parallelization scheme, by defining a `TATAMI_CUSTOM_PARALLEL` function-like macro. 
 *   This should accept the `fun`, `tasks` and `threads` arguments as below.
 * - `<thread>`, otherwise.
 *
 * @tparam parallel_ Whether the tasks should be run in parallel.
 * If `false`, no parallelization is performed and all tasks are run on the current thread.
 * @tparam Function_ Function to be applied for a contiguous range of tasks.
 * This should accept three arguments:
 * - `thread`, the thread number executing this task range.
 *   This will be passed as a `size_t`.
 * - `task_start`, the start index of the task range.
 *   This will be passed as a `Index_`.
 * - `task_length`, the number of tasks in the task range.
 *   This will be passed as a `Index_`.
 * @tparam Index_ Integer type for the number of tasks.
 *
 * @param fun Function that executes a contiguous range of tasks.
 * @param tasks Number of tasks.
 * @param threads Number of threads.
 */
template<bool parallel_ = true, class Function_, typename Index_>
void parallelize(Function_ fun, Index_ tasks, size_t threads) {
    if constexpr(parallel_) {
        if (threads > 1) {
#ifndef TATAMI_CUSTOM_PARALLEL
            Index_ worker_size = (tasks / threads) + (tasks % threads > 0); // Ceiling of an integer division.
            threads = (tasks / worker_size) + (tasks % worker_size > 0); // Set the actual number of required threads.
            std::vector<std::string> errors(threads);

#if defined(_OPENMP)
            #pragma omp parallel for num_threads(threads)
            for (size_t t = 0; t < threads; ++t) {
                Index_ start = worker_size * t; // Will not overflow due to the above recomputation of 'threads'.
                Index_ remaining = tasks - start; // Must be positive, as otherwise 'tasks % worker_size = 0' and the iteration wouldn't even get here.

                try {
                    fun(t, start, std::min(remaining, worker_size)); // Use 'remaining' to avoid potential overflow from computing 'end = start + worker_size'.
                } catch (std::exception& e) {
                    errors[t] = e.what();
                } catch (...) {
                    errors[t] = "unknown error in thread " + std::to_string(t);
                }
            }

#else
            Index_ first = 0;
            std::vector<std::thread> workers;
            workers.reserve(threads);

            for (size_t t = 0; t < threads && first < tasks; ++t) {
                Index_ remaining = tasks - first;
                Index_ len = std::min(remaining, worker_size);
                workers.emplace_back([&fun,&errors](int t, Index_ first, Index_ len) -> void {
                    try {
                        fun(t, first, len);
                    } catch (std::exception& e) {
                        errors[t] = e.what();
                    } catch (...) {
                        errors[t] = "unknown error in thread " + std::to_string(t);
                    }
                }, t, first, len);
                first += len;
            }

            for (auto& wrk : workers) {
                wrk.join();
            }
#endif

            for (const auto& e : errors) {
                if (!e.empty()) {
                    throw std::runtime_error(e);
                }
            }

#else
            TATAMI_CUSTOM_PARALLEL(std::move(fun), tasks, threads);
#endif
            return;
        }
    }

    fun(0, 0, tasks);
    return;
}

}

#endif
