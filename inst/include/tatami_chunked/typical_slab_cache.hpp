#ifndef TATAMI_CHUNKED_TYPICAL_SLAB_CACHE_HPP
#define TATAMI_CHUNKED_TYPICAL_SLAB_CACHE_HPP

#include <memory>
#include <type_traits>

#include "OracleSlabCache.hpp"
#include "SubsettedOracleSlabCache.hpp"
#include "LruSlabCache.hpp"

/**
 * @file typical_slab_cache.hpp
 * @brief Classes for typical caching of slabs. 
 */

namespace tatami_chunked {

/**
 * @brief Typical options for cached slab extraction.
 *
 * These options are usually relevant to any matrix representation that stores its data in a regular grid of rectangular chunks.
 * We define a "slab" as the set of chunks required to read a row/column (or a contiguous block/indexed subset thereof) during a `tatami::DenseExtractor::fetch()` or `tatami::SparseExtractor::Fetch()` call.
 * We aim to cache one or more complete slabs; this means that we can re-use the cached chunks when adjacent rows/columns are requested, rather than re-reading them (e.g., from disk).
 */
struct TypicalSlabCacheOptions {
    /**
     * Size of the in-memory cache in bytes.
     * Larger caches improve access speed at the cost of memory usage.
     * Small values may be ignored if `require_minimum_cache` is `true`.
     */
    size_t maximum_cache_size = 100000000;

    /**
     * Whether to automatically enforce a minimum size for the cache, regardless of `maximum_cache_size`.
     * This minimum is chosen to ensure that a single slab can be retained in memory,
     * so that the same chunks are not repeatedly re-read when iterating over consecutive rows/columns of the matrix.
     */
    bool require_minimum_cache = true;
};

/**
 * @brief Workspace for typical slab extraction.
 *
 * Implements a workspace to initialize the slab caches (i.e., the `LruSlabCache` and `OracleSlabCache`) for extraction from a chunked matrix representation.
 * This is intended to be a member of an `Extractor` class, allowing extraction to switch between different caches, e.g., when `ExtractorBase::set_oracle()` is called.
 * It also handles the calculation of various cache size statistics.
 *
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Slab_ Class that contains a single slab.
 * @tparam subset_ Whether to report the requested subset from each slab.
 * Only relevant when an oracle is available.
 */
template<typename Index_, class Slab_, bool subset_ = false>
struct TypicalSlabCacheWorkspace {
    /**
     * Default constructor.
     */
    TypicalSlabCacheWorkspace() = default;

    /**
     * @param primary_length Length of the primary dimension of each slab.
     * The primary dimension contains the elements to be extracted from each cached slab.
     * For example, if we were iterating through rows of a matrix, `primary_length` would be the number of rows spanned by each slab.
     * @param secondary_length Length of the secondary dimension of each slab.
     * This is the dimension that is not the primary, e.g., if we were iterating through matrix rows, the `secondary_length` would be the number of columns spanned by each slab.
     * @param cache_size_in_elements Total size of the cache in terms of the number of elements.
     * This is usually derived from `TypicalSlabCacheOptions::maximum_cache_size` and the size of each element.
     * @param require_minimum_cache Whether to enforce a minimum size of the cache for efficient extraction, see `TypicalSlabCacheOptions` for details.
     */
    TypicalSlabCacheWorkspace(Index_ primary_length, Index_ secondary_length, size_t cache_size_in_elements, bool require_minimum_cache) : primary_length(primary_length) {
        slab_size_in_elements = static_cast<size_t>(primary_length) * static_cast<size_t>(secondary_length);
        num_slabs_in_cache = (slab_size_in_elements ? cache_size_in_elements / slab_size_in_elements : 1);

        if (num_slabs_in_cache == 0 && require_minimum_cache) {
            num_slabs_in_cache = 1;
        }

        if (num_slabs_in_cache) {
            lru_cache.reset(new LruSlabCache<Index_, Slab_>(num_slabs_in_cache));
        }
    }

public:
    /**
     * Length of the primary dimension of each slab.
     */
    Index_ primary_length;

    /**
     * Size of each slab, in terms of the number of elements.
     */
    size_t slab_size_in_elements;

    /**
     * Number of slabs that can fit in the cache.
     */
    size_t num_slabs_in_cache;

    /**
     * Cache of least recently used slabs.
     * This is only allocated if `num_slabs_in_cache` is positive and `oracle_cache` is NULL.
     */
    std::unique_ptr<LruSlabCache<Index_, Slab_> > lru_cache;

    /**
     * Type of the oracle slab cache, depending on whether `subset_ = true`.
     */
    typedef typename std::conditional<subset_, SubsettedOracleSlabCache<Index_, Index_, Slab_>, OracleSlabCache<Index_, Index_, Slab_> >::type OracleCache;

    /**
     * Cache of to-be-used slabs, based on an `Oracle`'s predictions.
     * This may be NULL, see `set_oracle()` for more details.
     */
    std::unique_ptr<OracleCache> oracle_cache;

public:
    /**
     * Set up the oracle cache.
     * This will only have an effect if the number of slabs in the cache is greater than 1,
     * otherwise the predictions have no effect on the choice of slab to retain.
     * Callers should check for a non-NULL `oracle_cache` before attempting to use it,
     * otherwise they should use the `lru_cache` (provided `num_slabs_in_cache > 0`).
     *
     * @param o Oracle to provide predictions for subsequent accesses on the primary dimension.
     */
    void set_oracle(std::unique_ptr<tatami::Oracle<Index_> > o) {
        // The oracle won't have any effect if fewer than one slab can be cached.
        if (num_slabs_in_cache > 1) {
            size_t max_predictions = static_cast<size_t>(num_slabs_in_cache) * primary_length * 2; // double the cache size, basically.
            oracle_cache.reset(new OracleCache(std::move(o), max_predictions, num_slabs_in_cache));
            lru_cache.reset();
        }
    }
};

}

#endif
