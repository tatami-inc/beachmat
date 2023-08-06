#ifndef TATAMI_CHUNKED_LRU_SLAB_CACHE_HPP
#define TATAMI_CHUNKED_LRU_SLAB_CACHE_HPP

#include <unordered_map>
#include <list>

/**
 * @file LruSlabCache.hpp
 * @brief Create a LRU cache of slabs.
 */

namespace tatami_chunked {

/**
 * @tparam Id_ Type of cache identifier, typically integer.
 * @tparam Slab_ Class for a single slab.
 *
 * @brief Least-recently-used cache for slabs.
 *
 * Implements a least-recently-used (LRU) cache, typically containing one or more "slabs" from a chunked matrix representation.
 * Each slab is defined as the set of chunks required to read a row/column (or a contiguous block/indexed subset thereof) during iteration through a `tatami::Matrix`.
 * The LRU cache can be used for chunked `tatami::Matrix` representations where the data is costly to load (e.g., from file) and no oracle is provided to predict future accesses.
 * In such cases, chunks of data can be loaded and cached such that any possible future request for an already-loaded slab will just fetch it from cache.
 */
template<typename Id_, class Slab_> 
class LruSlabCache {
private:
    typedef std::pair<Slab_, Id_> Element;
    std::list<Element> cache_data;
    std::unordered_map<Id_, typename std::list<Element>::iterator> cache_exists;
    size_t max_slabs; 

public:
    /**
     * @param m Maximum number of slabs to store.
     */
    LruSlabCache(size_t m) : max_slabs(m) {}

public:
    /**
     * @tparam Cfunction_ Function to create a new `Slab_` object.
     * @tparam Pfunction_ Function to populate a `Slab_` object with the contents of a slab.
     *
     * @param id Identifier for the cached slab.
     * This is typically defined as the index of the slab on the iteration dimension.
     * For example, if each chunk takes up 10 rows, attempting to access row 21 would require retrieval of slab 2.
     * @param create Function that accepts no arguments and returns a `Slab_` object.
     * @param populate Function that accepts a slab ID and a reference to a `Slab_` object,
     * and populates the latter with the contents of the former.
     * 
     * @return Reference to a slab.
     * If the slab already exists in the cache, it is returned directly.
     * If the slab does not exist and there is still space in the cache, a new slab is created and populated with the contents of slab `id`.
     * If the slab does not exist and there is no space in the cache, the least recently used slab is evicted and its `Slab_` is populated with the contents of slab `id`.
     */
    template<class Cfunction_, class Pfunction_>
    const Slab_& find(Id_ id, Cfunction_ create, Pfunction_ populate) {
        if (max_slabs == 1) {
            // Minor optimization if there's just one slab, in which case we can
            // skip the search and the list splice if there's a hit.
            if (!cache_data.empty()) {
                const auto& solo = cache_data.front();
                if (solo.second == id) {
                    return solo.first;
                }
            }
        } else {
            auto it = cache_exists.find(id);
            if (it != cache_exists.end()) {
                auto chosen = it->second;
                cache_data.splice(cache_data.end(), cache_data, chosen); // move to end.
                return chosen->first;
            } 
        }

        typename std::list<Element>::iterator location;
        if (cache_data.size() < max_slabs) {
            cache_data.emplace_back(create(), id);
            location = std::prev(cache_data.end());
        } else {
            location = cache_data.begin();
            cache_exists.erase(location->second);
            location->second = id;
            cache_data.splice(cache_data.end(), cache_data, location); // move to end.
        }
        cache_exists[id] = location;

        auto& slab = location->first;
        populate(id, slab);
        return slab;
    }
};

}

#endif
