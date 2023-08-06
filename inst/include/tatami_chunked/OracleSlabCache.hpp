#ifndef TATAMI_CHUNKED_ORACLE_SLAB_CACHE_HPP
#define TATAMI_CHUNKED_ORACLE_SLAB_CACHE_HPP

#include <unordered_map>
#include <vector>
#include <list>
#include "tatami/tatami.hpp"

/**
 * @file OracleSlabCache.hpp
 * @brief Create a oracle-aware cache for slabs.
 */

namespace tatami_chunked {

/**
 * @brief Oracle-aware cache for slabs.
 *
 * @tparam Id_ Type of slab identifier, typically integer.
 * @tparam Index_ Type of row/column index produced by the oracle.
 * @tparam Slab_ Class for a single slab.
 *
 * Implement an oracle-aware cache for slabs.
 * Each slab is defined as the set of chunks required to read a row/column (or a contiguous block/indexed subset thereof) during iteration through a `tatami::Matrix`.
 * This cache can be used for `Matrix` representations where the data is costly to load (e.g., from file) and a `tatami::Oracle` is provided to predict future accesses.
 * In such cases, chunks of data can be loaded and cached such that any possible future request for an already-loaded slab will just fetch it from cache.
 */
template<typename Id_, typename Index_, class Slab_> 
class OracleSlabCache {
    tatami::OracleStream<Index_> prediction_stream;
    size_t max_predictions;
    size_t max_slabs;

private:
    std::list<Slab_> slab_cache, tmp_cache, free_cache;

    typedef typename std::list<Slab_>::iterator cache_iterator;
    std::unordered_map<Id_, std::pair<Index_, cache_iterator> > slab_exists, past_exists;

    std::vector<std::pair<Index_, Index_> > predictions_made;
    size_t predictions_fulfilled = 0;

    std::vector<Slab_*> slab_pointers;

    std::vector<std::pair<Index_, cache_iterator*> > unassigned_slabs;
    std::vector<std::pair<Id_, Index_> > slabs_to_populate; 

public:
    /**
     * @param oracle Pointer to an `tatami::Oracle` to be used for predictions.
     * @param per_iteration Maximum number of predictions to make per iteration.
     * @param num_slabs Maximum number of slabs to store.
     */
    OracleSlabCache(std::unique_ptr<tatami::Oracle<Index_> > oracle, size_t per_iteration, size_t num_slabs) :
        prediction_stream(std::move(oracle)), 
        max_predictions(per_iteration),
        max_slabs(num_slabs)
    {
        slab_exists.reserve(max_slabs);
        past_exists.reserve(max_slabs);

        predictions_made.reserve(max_predictions);

        slab_pointers.reserve(max_slabs);
        unassigned_slabs.reserve(max_slabs);
        slabs_to_populate.reserve(max_slabs);
    } 

    /**
     * @cond
     */
    // For testing only.
    OracleSlabCache() = default;
    /**
     * @endcond
     */

private:
    std::pair<const Slab_*, Index_> fetch(size_t i) const {
        const auto& current = predictions_made[i];
        return std::pair<const Slab_*, Index_>(slab_pointers[current.first], current.second);
    }

public:
    /**
     * Fetch the next slab according to the stream of predictions provided by the `tatami::Oracle`.
     *
     * @tparam Ifunction_ Function to identify the slab containing each predicted row/column.
     * @tparam Cfunction_ Function to create a new slab.
     * @tparam Pfunction_ Function to populate zero, one or more slabs with their contents.
     *
     * @param identify Function that accepts an `i`, an `Index_` containing the predicted row/column index.
     * This should return a pair containing (1) the identifier of the slab containing `i`, and (2) the index of row/column `i` inside that slab.
     * This is typically defined as the index of the slab on the iteration dimension.
     * For example, if each chunk takes up 10 rows, attempting to access row 21 would require retrieval of slab 2 and an offset of 1.
     * @param create Function that accepts no arguments and returns a `Slab_` object with sufficient memory to hold a slab's contents when used in `populate()`.
     * This may also return a default-constructed `Slab_` object if the allocation is done dynamically per slab in `populate()`.
     * @param populate Function that accepts two arguments, `slabs_in_need` and `slab_data`.
     * (1) `slabs_in_need` is a `const std::vector<std::pair<Id_, Index_> >&` specifying the slabs to be populated.
     * The first `Id_` element of each pair contains the slab identifier, i.e., the first element returned by the `identify` function.
     * The second `Index_` element specifies the index in `slab_data` in which to store the contents of each slab.
     * (2) `slab_data` is a `std::vector<Slab_*>&` containing pointers to the cached slab contents to be populated.
     * This function should iterate over the `slabs_in_need` and populate the corresponding entries in `slab_data`.
     *
     * @return Pair containing (1) a pointer to a slab's contents and (2) the index of the next predicted row/column inside the retrieved slab.
     */
    template<class Ifunction_, class Cfunction_, class Pfunction_>
    std::pair<const Slab_*, Index_> next(Ifunction_ identify, Cfunction_ create, Pfunction_ populate) {
        if (predictions_made.size() > predictions_fulfilled) {
            return fetch(predictions_fulfilled++);
        }

        predictions_made.clear();
        size_t used = 0;

        // Iterators in the unordered_map should remain valid after swapping the containers, 
        // see https://stackoverflow.com/questions/4124989/does-stdvectorswap-invalidate-iterators
        tmp_cache.swap(slab_cache);

        past_exists.swap(slab_exists);
        slab_exists.clear();

        slab_pointers.clear();
        unassigned_slabs.clear();
        slabs_to_populate.clear();

        for (size_t p = 0; p < max_predictions; ++p) {
            Index_ current;
            if (!prediction_stream.next(current)) {
                break;
            }

            auto slab_id = identify(current);
            auto curslab = slab_id.first;
            auto curindex = slab_id.second;

            auto it = slab_exists.find(curslab);
            if (it != slab_exists.end()) {
                predictions_made.emplace_back((it->second).first, curindex);

            } else if (used < max_slabs) {
                auto past = past_exists.find(curslab);
                if (past != past_exists.end()) {
                    auto sIt = (past->second).second;
                    slab_cache.splice(slab_cache.end(), tmp_cache, sIt);
                    slab_pointers.push_back(&(*sIt));
                    slab_exists[curslab] = std::make_pair(used, sIt);
                    
                } else {
                    if (free_cache.empty()) {
                        // We might be able to recycle an existing slab from tmp_cache 
                        // to populate 'curslab'... but we don't know if we can do so at
                        // this moment, as those slabs might be needed by later predictions.
                        // So we just defer the creation of a new slab until we've run 
                        // through the set of predictions for this round.
                        auto ins = slab_exists.insert(std::make_pair(curslab, std::make_pair(used, slab_cache.end())));
                        unassigned_slabs.emplace_back(used, &(ins.first->second.second));
                        slab_pointers.push_back(NULL);

                    } else {
                        auto sIt = free_cache.begin();
                        slab_cache.splice(slab_cache.end(), free_cache, sIt);
                        slab_pointers.push_back(&(*sIt));
                        slab_exists[curslab] = std::make_pair(used, sIt);
                    }

                    slabs_to_populate.emplace_back(curslab, used);
                }

                predictions_made.emplace_back(used, curindex);
                ++used;

            } else {
                prediction_stream.back();
                break;
            }
        }

        while (!unassigned_slabs.empty()) {
            cache_iterator it;
            if (!tmp_cache.empty()) {
                it = tmp_cache.begin();
                slab_cache.splice(slab_cache.end(), tmp_cache, it);
            } else {
                slab_cache.emplace_back(create());
                it = slab_cache.end();
                --it;
            }

            auto& last = unassigned_slabs.back();
            slab_pointers[last.first] = &(*it);

            // This changes the value in the slab_exists map without having to do a look-up, see:
            // https://stackoverflow.com/questions/16781886/can-we-store-unordered-maptiterator
            *(last.second) = it; 

            unassigned_slabs.pop_back();
        }

        while (!tmp_cache.empty()) {
            free_cache.splice(free_cache.end(), tmp_cache, tmp_cache.begin());
        }

        if (!slabs_to_populate.empty()) {
            populate(slabs_to_populate, slab_pointers);
        }

        // Well, because we just used one.
        predictions_fulfilled = 1;
        return fetch(0);
    }
};

}

#endif
