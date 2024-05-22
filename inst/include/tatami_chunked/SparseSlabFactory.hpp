#ifndef TATAMI_CHUNKED_SPARSE_SLAB_FACTORY_HPP
#define TATAMI_CHUNKED_SPARSE_SLAB_FACTORY_HPP

#include <vector>
#include "SlabCacheStats.hpp"

/**
 * @file SparseSlabFactory.hpp 
 * @brief Factory for sparse slabs.
 */

namespace tatami_chunked {

/**
 * @brief Factory for sparse slabs.
 *
 * An instance of this class allocates two memory pools (for values and indices) during its construction.
 * Each slab simply contains pointers into this pool at regular offsets.
 * The idea is to allocate large blocks for caching in, e.g., `LruSlabCache`, improving performance by avoiding repeated requests to the allocator.
 * This also reduces fragmentation that could increase memory usage beyond the expected cache size.
 *
 * @tparam Value_ Type of the data in each slab.
 * @tparam Index_ Integer type of the indices in each slab.
 * @tparam Count_ Integer type for counting structural non-zeros.
 * This should be large enough to store the extent of the non-target dimension of the slab.
 */
template<typename Value_, typename Index_, typename Count_ = Index_>
struct SparseSlabFactory {
    /**
     * @param target_dim Extent of the target dimension of the slab,
     * i.e., the dimension that is indexed into.
     * @param non_target_dim Extent of the non-target dimension of the slab.
     * @param slab_size Size of the slab.
     * This should be equal to the product of `target_dim` and `non_target_dim`.
     * @param max_slabs Maximum number of slabs.
     * @param needs_value Whether the values of the structural non-zeros should be cached.
     * @param needs_index Whether the indices of the structural non-zeros should be cached.
     */
    SparseSlabFactory(size_t target_dim, size_t non_target_dim, size_t slab_size, size_t max_slabs, bool needs_value, bool needs_index) : 
        my_target_dim(target_dim),
        my_non_target_dim(non_target_dim),
        my_slab_size(slab_size),
        my_needs_value(needs_value),
        my_needs_index(needs_index),
        my_number_pool(max_slabs * target_dim)
    {
        size_t total_size = max_slabs * slab_size;
        if (needs_value) {
            my_value_pool.resize(total_size);
        }
        if (needs_index) {
            my_index_pool.resize(total_size);
        }
    }

    /**
     * Overload that computes `slab_size` automatically.
     * @param target_dim Extent of the target dimension of the slab.
     * @param non_target_dim Extent of the non-target dimension of the slab.
     * @param max_slabs Maximum number of slabs.
     * @param needs_value Whether the values of the structural non-zeros should be cached.
     * @param needs_index Whether the indices of the structural non-zeros should be cached.
     */
    SparseSlabFactory(size_t target_dim, size_t non_target_dim, size_t max_slabs, bool needs_value, bool needs_index) : 
        SparseSlabFactory(target_dim, non_target_dim, target_dim * non_target_dim, max_slabs, needs_value, needs_index) {}

    /**
     * Overload that takes the relevant statistics from a `SlabCacheStats` object.
     * @param target_dim Extent of the target dimension of the slab.
     * @param non_target_dim Extent of the non-target dimension of the slab.
     * @param stats Slab statistics, computed from `target_dim` and `non_target_dim`.
     * @param needs_value Whether the values of the structural non-zeros should be cached.
     * @param needs_index Whether the indices of the structural non-zeros should be cached.
     */
    SparseSlabFactory(size_t target_dim, size_t non_target_dim, const SlabCacheStats& stats, bool needs_value, bool needs_index) : 
        SparseSlabFactory(target_dim, non_target_dim, stats.slab_size_in_elements, stats.max_slabs_in_cache, needs_value, needs_index) {}

    /**
     * @cond
     */
    // Delete the copy constructors as we're passing out pointers.
    SparseSlabFactory(const SparseSlabFactory&) = delete;
    SparseSlabFactory& operator=(const SparseSlabFactory&) = delete;

    // Move constructors are okay though.
    SparseSlabFactory(SparseSlabFactory&&) = default;
    SparseSlabFactory& operator=(SparseSlabFactory&&) = default;
    /**
     * @endcond
     */

private:
    size_t my_offset_slab = 0, my_offset_number = 0;
    size_t my_target_dim, my_non_target_dim, my_slab_size;
    bool my_needs_value, my_needs_index;
    std::vector<Value_> my_value_pool;
    std::vector<Index_> my_index_pool;
    std::vector<Count_> my_number_pool;

public:
    /**
     * @brief Sparse slab.
     */
    struct Slab {
        /**
         * Vector of pointers of length equal to `target_dim`. 
         * Each pointer corresponds to an element of the target dimension of the slab,
         * and refers to an array with `non_target_dim` addressable elements.
         * Each pointer should be used to store the values of the structural non-zeros for that dimension element.
         *
         * Alternatively, this vector may be empty if `needs_value = false` in the `SparseSlabFactory` constructor.
         */
        std::vector<Value_*> values;

        /**
         * Vector of pointers of length equal to `target_dim`. 
         * Each pointer corresponds to an element of the target dimension of the slab,
         * and refers to an array with `non_target_dim` addressable elements.
         * Each pointer should be used to store the indices of the structural non-zeros for that dimension element.
         *
         * Alternatively, this vector may be empty if `needs_index = false` in the `SparseSlabFactory` constructor.
         */
        std::vector<Index_*> indices;

        /**
         * Pointer to an array with `target_dim` addressable elements.
         * Each value stores the number of non-zero elements for each element of the target dimension of the slab,
         * i.e., the number of entries that are filled in the corresponding arrays of `values` and `indices`.
         * On creation, all entries of this array are set to zero.
         */
        Count_* number = NULL;
    };

    /**
     * Create a new slab, i.e., designate a portion of each memory pool for use through the returned pointers.
     * This should not be called more than `max_slabs` times.
     *
     * @return Slab containing pointers to the relevant memory pools.
     */
    Slab create() {
        Slab output;
        output.number = my_number_pool.data() + my_offset_number;
        my_offset_number += my_target_dim;

        if (my_needs_value) {
            output.values.reserve(my_target_dim);
            auto vptr = my_value_pool.data() + my_offset_slab;
            for (size_t p = 0; p < my_target_dim; ++p, vptr += my_non_target_dim) {
                output.values.push_back(vptr);
            }
        }

        if (my_needs_index) {
            output.indices.reserve(my_target_dim);
            auto iptr = my_index_pool.data() + my_offset_slab;
            for (size_t p = 0; p < my_target_dim; ++p, iptr += my_non_target_dim) {
                output.indices.push_back(iptr);
            }
        }

        my_offset_slab += my_slab_size;
        return output;
    }
};

}

#endif
