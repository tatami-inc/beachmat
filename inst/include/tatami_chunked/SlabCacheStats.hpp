#ifndef TATAMI_CHUNKED_SLAB_CACHE_STATS_HPP
#define TATAMI_CHUNKED_SLAB_CACHE_STATS_HPP

#include <algorithm>

/**
 * @file SlabCacheStats.hpp
 * @brief Slab cache statistics.
 */

namespace tatami_chunked {

/**
 * @brief Statistics for slab caching.
 *
 * This computes the slab size and the number of slabs to be cached, given the dimensions of the slab and the cache size in bytes.
 * The assumption is that all slabs are of the same shape, partitioning the matrix into regular intervals along the target dimension.
 * Developers should check out `CustomDenseChunkedMatrix` for some usage examples.
 */
struct SlabCacheStats {
    /**
     * Size of each slab, in terms of the number of data elements.
     */
    size_t slab_size_in_elements;

    /**
     * Number of slabs that can fit in the cache.
     * This is used as `max_slabs` in `LruSlabCache`, `OracularSlabCache` and friends.
     */
    size_t max_slabs_in_cache;

    /**
     * @param target_length Length of the target dimension of each slab.
     * For example, if we were iterating through rows of a matrix, `target_length` would be the number of rows spanned by each slab.
     * @param non_target_length Length of the non-target dimension of each slab.
     * For example, if we were iterating through matrix rows, the `non_target_length` would be the number of columns spanned by each slab.
     * @param target_num_slabs Number of slabs required to span the full extent of the target dimension of the matrix.
     * This is used as an upper bound on the value of `max_slabs_in_cache`.
     * For example, if we were iterating through matrix rows, the `target_num_slabs` would be the number of slabs required to span all rows of the matrix.
     * @param cache_size_in_elements Total size of the cache, in terms of the number of data elements.
     * @param require_minimum_cache Whether to enforce a minimum size of the cache for efficient extraction of consecutive dimension elements, even if it exceeds `cache_size_in_elements`.
     */
    SlabCacheStats(size_t target_length, size_t non_target_length, size_t target_num_slabs, size_t cache_size_in_elements, bool require_minimum_cache) :
        slab_size_in_elements(target_length * non_target_length),
        max_slabs_in_cache(compute_max_slabs_in_cache(slab_size_in_elements, target_num_slabs, cache_size_in_elements, require_minimum_cache))
    {}

    /**
     * @param target_length Length of the target dimension of each slab.
     * For example, if we were iterating through rows of a matrix, `target_length` would be the number of rows spanned by each slab.
     * @param non_target_length Length of the non-target dimension of each slab.
     * For example, if we were iterating through matrix rows, the `non_target_length` would be the number of columns spanned by each slab.
     * @param target_num_slabs Number of slabs required to span the full extent of the target dimension of the matrix.
     * This is used as an upper bound on the value of `max_slabs_in_cache`.
     * For example, if we were iterating through matrix rows, the `target_num_slabs` would be the number of slabs required to span all rows of the matrix.
     * @param cache_size_in_bytes Total size of the cache, in terms of the number of bytes.
     * @param element_size Size of each data element in the cache, in bytes.
     * This may be zero, e.g., when neither the value nor the index are required during sparse extraction.
     * @param require_minimum_cache Whether to enforce a minimum size of the cache for efficient extraction of consecutive dimension elements, even if it exceeds `cache_size_in_bytes`.
     */
    SlabCacheStats(size_t target_length, size_t non_target_length, size_t target_num_slabs, size_t cache_size_in_bytes, size_t element_size, bool require_minimum_cache) :
        slab_size_in_elements(target_length * non_target_length),
        max_slabs_in_cache([&]() {
            if (element_size == 0) {
                return target_num_slabs;
            } else {
                return compute_max_slabs_in_cache(slab_size_in_elements, target_num_slabs, cache_size_in_bytes / element_size, require_minimum_cache); 
            }
        }())
    {}

private:
    static size_t compute_max_slabs_in_cache(size_t slab_size_in_elements, size_t num_slabs, size_t cache_size_in_elements, bool require_minimum_cache) {
        if (slab_size_in_elements == 0) {
            return num_slabs;
        }

        auto tmp = cache_size_in_elements / slab_size_in_elements;
        if (tmp == 0 && require_minimum_cache) {
            return 1;
        } 

        return std::min(tmp, num_slabs);
    }
};

}

#endif
