#ifndef TATAMI_CHUNKED_DENSE_SLAB_FACTORY_HPP
#define TATAMI_CHUNKED_DENSE_SLAB_FACTORY_HPP

#include <vector>
#include "SlabCacheStats.hpp"

/**
 * @file DenseSlabFactory.hpp 
 * @brief Factory for dense slabs.
 */

namespace tatami_chunked {

/**
 * @brief Factory for dense slabs.
 *
 * An instance of this class allocates a single memory pool during its construction.
 * Each slab is simply a pointer into this pool at regular offsets.
 * The idea is to allocate one large block for caching in, e.g., `LruSlabCache`, improving performance by avoiding repeated requests to the allocator.
 * This also reduces fragmentation that could increase memory usage beyond the expected cache size.
 *
 * @tparam Value_ Type of the data in each slab.
 */
template<typename Value_>
struct DenseSlabFactory {
    /**
     * @param slab_size Size of the slab, in terms of data elements.
     * @param max_slabs Maximum number of slabs.
     */
    DenseSlabFactory(size_t slab_size, size_t max_slabs) : my_slab_size(slab_size), my_pool(max_slabs * slab_size) {}

    /**
     * @param stats Slab cache statistics.
     */
    DenseSlabFactory(const SlabCacheStats& stats) : DenseSlabFactory(stats.slab_size_in_elements, stats.max_slabs_in_cache) {}

    /**
     * @cond
     */
    // Delete the copy constructors as we're passing out pointers.
    DenseSlabFactory(const DenseSlabFactory&) = delete;
    DenseSlabFactory& operator=(const DenseSlabFactory&) = delete;

    // Move constructors are okay though.
    DenseSlabFactory(DenseSlabFactory&&) = default;
    DenseSlabFactory& operator=(DenseSlabFactory&&) = default;
    /**
     * @endcond
     */

private:
    size_t my_offset = 0, my_slab_size;
    std::vector<Value_> my_pool;

public:
    /**
     * @brief Dense slab.
     */
    struct Slab {
        /**
         * Pointer to a buffer with `slab_size` addressable data elements (as used in the `DenseSlabFactory` constructor).
         */
        Value_* data = NULL;
    };

    /**
     * Create a new slab, i.e., designate a portion of the memory pool for use through the returned pointer.
     * This should not be called more than `max_slabs` times.
     *
     * @return Slab containing a pointer to the memory pool.
     */
    Slab create() {
        Slab output;
        output.data = my_pool.data() + my_offset;
        my_offset += my_slab_size;
        return output;
    }
};

}

#endif
