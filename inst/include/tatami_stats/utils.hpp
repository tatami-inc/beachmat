#ifndef TATAMI_STATS__UTILS_HPP
#define TATAMI_STATS__UTILS_HPP

#include <vector>
#include <algorithm>

/**
 * @file utils.hpp
 *
 * @brief Utilities for computing matrix statistics.
 */

namespace tatami_stats {

/**
 * Count the total number of groups, typically for per-group memory allocations.
 *
 * @tparam Group_ Integer type for the group assignments.
 * @tparam Index_ Integer type for the number of observations.
 *
 * @param[in] group Pointer to an array of group assignments per observation.
 * Each assignment should be an integer in `[0, G)` where `G` is the total number of groups.
 * @param n Number of observations, i.e., the length of the array referenced by `group`.
 *
 * @return Total number of groups, i.e., `G`.
 * Note that not all groups may actually have non-zero occurrences in `group`.
 */
template<typename Group_, typename Size_>
size_t total_groups(const Group_* group, Size_ n) {
    if (n) {
        return static_cast<size_t>(*std::max_element(group, group + n)) + 1;
    } else {
        return 0;
    }
}

/**
 * Count the occurrences of each group.
 *
 * @tparam Group_ Integer type for the group assignments.
 * @tparam Index_ Integer type for the number of observations.
 *
 * @param[in] group Pointer to an array of group assignments per observation.
 * Each assignment should be an integer in `[0, G)` where `G` is the total number of groups.
 * @param n Number of observations, i.e., the length of the array referenced by `group`.
 *
 * @return Vector of length equal to `G`, containing the number of occurrences of each group.
 */
template<typename Group_, typename Size_>
std::vector<Size_> tabulate_groups(const Group_* group, Size_ n) {
    auto ngroups = total_groups(group, n);
    std::vector<Size_> group_sizes(ngroups);
    for (Size_ r = 0; r < n; ++r) {
        ++(group_sizes[group[r]]);
    }
    return group_sizes;
}

/**
 * @brief Local output buffer for running calculations.
 *
 * A typical parallelization scenario involves dividing the set of objective vectors into contiguous blocks, where each thread operates on a block at a time.
 * However, in running calculations, an entire block's statistics are updated when its corresponding thread processes an observed vector.
 * If these statistics are stored in a global buffer, false sharing at the boundaries of the blocks can result in performance degradation. 
 *
 * To avoid this, the `LocalOutputBuffer` class provides thread-local storage for output statistics. 
 * Once the calculations are finished per thread, callers should use `transfer()` to transfer the local statistics to the global buffer.
 * The exception is that of the first thread, which is allowed to directly write to the global output buffer.
 *
 * @tparam Output_ Type of the result.
 */
template<typename Output_>
class LocalOutputBuffer {
public:
    /**
     * @tparam Index_ Type of the start index and length.
     * @param thread Identity of the thread, starting from zero to the total number of threads.
     * @param start Index of the first objective vector in the contiguous block for this thread.
     * @param length Number of objective vectors in the contiguous block for this thread.
     * @param[out] output Pointer to the global output buffer.
     */
    template<typename Index_>
    LocalOutputBuffer(size_t thread, Index_ start, Index_ length, Output_* output) : my_output(output + start), use_local(thread > 0), my_buffer(use_local ? length : 0) {
        if (!use_local) {
            // Setting to zero to match the initial behavior of 'my_buffer' when 'use_local = true'.
            std::fill_n(my_output, length, static_cast<Output_>(0));
        }
    }

    /**
     * Default constructor.
     */
    LocalOutputBuffer() = default;

    /**
     * @return Pointer to an output buffer to use for this thread.
     * This contains at least `length` addressable elements (see the argument of the same name in the constructor). 
     */
    Output_* data() {
        return (use_local ? my_buffer.data() : my_output);
    }

    /**
     * Transfer results from the local buffer to the global buffer (i.e., `output` in the constructor).
     */
    void transfer() {
        if (use_local) {
            std::copy(my_buffer.begin(), my_buffer.end(), my_output);
        }
    }

private:
    Output_* my_output = NULL;
    bool use_local = false;
    std::vector<Output_> my_buffer;
};

}

#endif
