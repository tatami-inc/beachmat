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
 * Each slab is defined as the set of chunks required to read an element of the target dimension (or a contiguous block/indexed subset thereof) from a `tatami::Matrix`.
 * The LRU cache can be used for chunked `tatami::Matrix` representations where the data is costly to load (e.g., from file) and no oracle is provided to predict future accesses on the target dimension.
 * In such cases, chunks of data can be loaded and cached such that any possible future request for an already-loaded slab will just fetch it from cache.
 */
template<typename Id_, class Slab_> 
class LruSlabCache {
private:
    typedef std::pair<Slab_, Id_> Element;
    std::list<Element> my_cache_data;
    std::unordered_map<Id_, typename std::list<Element>::iterator> my_cache_exists;
    size_t my_max_slabs;
    Id_ my_last_id = 0;
    Slab_* my_last_slab = NULL;

public:
    /**
     * @param max_slabs Maximum number of slabs to store in the cache.
     */
    LruSlabCache(size_t max_slabs) : my_max_slabs(max_slabs) {}

    /**
     * Deleted as the cache holds persistent iterators.
     */
    LruSlabCache(const LruSlabCache&) = delete;

    /**
     * Deleted as the cache holds persistent iterators.
     */
    LruSlabCache& operator=(const LruSlabCache&) = delete;

    /**
     * @cond
     */
    // Iterators are guaranteed to be valid after move, see Notes in
    // https://en.cppreference.com/w/cpp/container/list/list
    // https://en.cppreference.com/w/cpp/container/list/operator%3D
    LruSlabCache& operator=(LruSlabCache&&) = default; 
    LruSlabCache(LruSlabCache&&) = default;

    // Might as well define this.
    ~LruSlabCache() = default;
    /**
     * @endcond
     */

public:
    /**
     * This method should only be called if `m > 0` in the constructor.
     *
     * @tparam Cfunction_ Function to create a new `Slab_` object.
     * @tparam Pfunction_ Function to populate a `Slab_` object with the contents of a slab.
     *
     * @param id Identifier for the cached slab.
     * This is typically defined as the index of the slab on the target dimension.
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
        if (id == my_last_id && my_last_slab) {
            return *my_last_slab;
        }
        my_last_id = id;

        auto it = my_cache_exists.find(id);
        if (it != my_cache_exists.end()) {
            auto chosen = it->second;
            my_cache_data.splice(my_cache_data.end(), my_cache_data, chosen); // move to end.
            my_last_slab = &(chosen->first);
            return chosen->first;
        } 

        typename std::list<Element>::iterator location;
        if (my_cache_data.size() < my_max_slabs) {
            my_cache_data.emplace_back(create(), id);
            location = std::prev(my_cache_data.end());
        } else {
            location = my_cache_data.begin();
            my_cache_exists.erase(location->second);
            location->second = id;
            my_cache_data.splice(my_cache_data.end(), my_cache_data, location); // move to end.
        }
        my_cache_exists[id] = location;

        auto& slab = location->first;
        populate(id, slab);
        my_last_slab = &slab;
        return slab;
    }

public:
    /**
     * @return Maximum number of slabs in the cache.
     */
    size_t get_max_slabs() const {
        return my_max_slabs;
    }

    /**
     * @return Number of slabs currently in the cache.
     */
    size_t get_num_slabs() const {
        return my_cache_data.size();
    }
};

}

#endif
