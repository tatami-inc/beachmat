#ifndef TATAMI_R_PARALLELIZE_HPP
#define TATAMI_R_PARALLELIZE_HPP

/**
 * @cond
 */
#ifdef TATAMI_R_PARALLELIZE_UNKNOWN
/**
 * @endcond
 */

#include "manticore/manticore.hpp"
#include <thread>
#include <cmath>
#include <vector>
#include <string>
#include <stdexcept>
#include <algorithm>

/**
 * @file parallelize.hpp
 *
 * @brief Safely parallelize for unknown matrices.
 */

namespace tatami_r {

/**
 * Retrieve a global `manticore::Executor` object for all **tatami_r** uses.
 * This function is only available if `TATAMI_R_PARALLELIZE_UNKNOWN` is defined.
 *
 * @return Reference to a global `manticore::Executor`.
 */
inline manticore::Executor& executor() {
    // This should end up resolving to a single instance, even across dynamically linked libraries:
    // https://stackoverflow.com/questions/52851239/local-static-variable-linkage-in-a-template-class-static-member-function
    static manticore::Executor mexec;
    return mexec;
}

/**
 * @tparam Function_ Function to be executed.
 *
 * @param njobs Number of jobs to be executed.
 * @param fun Function to run in each thread.
 * This is a lambda that should accept three arguments:
 * - Integer containing the thread ID
 * - Integer specifying the index of the first job to be executed in a thread.
 * - Integer specifying the number of jobs to be executed in a thread.
 * @param nthreads Number of threads to parallelize over.
 *
 * The series of integers from 0 to `njobs - 1` is split into `nthreads` contiguous ranges.
 * Each range is used as input to `fun` within the corresponding thread.
 * It is assumed that the execution of any given job is independent of the next.
 *
 * This function is only available if `TATAMI_R_PARALLELIZE_UNKNOWN` is defined.
 */ 
template<class Function_>
void parallelize(Function_ fun, size_t njobs, size_t nthreads) {
    if (nthreads == 1 || njobs == 1) {
        fun(0, 0, njobs);
        return;
    }

    auto& mexec = executor();
    mexec.initialize(nthreads, "failed to execute R command");

    size_t jobs_per_worker = std::ceil(static_cast<double>(njobs) / nthreads);
    size_t start = 0;

    std::vector<std::thread> runners;
    runners.reserve(nthreads);
    std::vector<std::string> errors(nthreads);

    for (size_t w = 0; w < nthreads; ++w) {
        size_t end = std::min(njobs, start + jobs_per_worker);
        if (start >= end) {
            mexec.finish_thread(false);
            continue;
        }

        runners.emplace_back([&](size_t id, size_t s, size_t l) -> void {
            try {
                fun(id, s, l);
            } catch (std::exception& x) {
                // No throw here, we need to make sure we mark the
                // thread as being completed so that the main loop can quit.
                errors[id] = x.what();
            }
            mexec.finish_thread();
        }, w, start, end - start);

        start += jobs_per_worker;
    }

    mexec.listen();
    for (auto& x : runners) {
        x.join();
    }

    for (auto err : errors) {
        if (!err.empty()) {
            throw std::runtime_error(err);
        }
    }
}

}

/**
 * @cond
 */
#endif
/**
 * @endcond
 */

#endif
