#ifndef TATAMI_CHUNKED_SUBSETTED_ORACLE_SLAB_CACHE_HPP
#define TATAMI_CHUNKED_SUBSETTED_ORACLE_SLAB_CACHE_HPP

#include <unordered_map>
#include <vector>
#include <list>
#include "tatami/tatami.hpp"

/**
 * @file SubsettedOracleSlabCache.hpp
 * @brief Create a oracle-aware cache with subsets.
 */

namespace tatami_chunked {

/**
 * Type of subset selection.
 * Used to determine the subsets to extract in `SubsettedOracleSlabCache::SubsetDetails`.
 * 
 * - `FULL`: all rows/columns. 
 * - `BLOCK`: a contiguous block of rows/columns. 
 * - `INDEX`: an indexed subset of rows/columns.
 */
enum class SubsetSelection : char { FULL, BLOCK, INDEX };

/**
 * @brief Oracle-aware cache for slabs, plus subsets.
 *
 * @tparam Id_ Type of slab identifier, typically integer.
 * @tparam Index_ Type of row/column index produced by the oracle.
 * @tparam Slab_ Class for a single slab.
 *
 * Implement an oracle-aware cache for slab subsets.
 * This is similar to the `OracleSlabCache` except that it remembers the subset of rows/columns that were requested for each slab.
 * Slab extractors can use this information to optimize the slab population process by ignoring rows/columns that are not needed.
 */
template<typename Id_, typename Index_, class Slab_> 
class SubsettedOracleSlabCache {
    tatami::OracleStream<Index_> prediction_stream;
    size_t max_predictions;
    size_t max_slabs;

public:
    /**
     * @brief Details on the subset to extract.
     */
    struct SubsetDetails {
        /**
         * Type of subset selection to extract from the slab.
         */
        SubsetSelection selection;

        /**
         * Start of the block within the slab to be extracted.
         * Only used if `selection` is set to `SubsetSelection::BLOCK`.
         */
        Index_ block_start;

        /**
         * Length of the block within the slab to be extracted.
         * Only used if `selection` is set to `SubsetSelection::BLOCK`.
         */
        Index_ block_length;

        /**
         * Indices within the slab to be extracted.
         * Guaranteed to be sorted and unique.
         * Only used if `selection` is set to `SubsetSelection::INDEX`.
         */
        std::vector<Index_> indices;

        /**
         * Mapping of indices to be extracted to their positions inside `indices`.
         * All values of `indices` are present as keys here where `mapping[indices[i]] = i`.
         * Only used if `selection` is set to `SubsetSelection::INDEX`.
         */
        std::unordered_map<Index_, Index_> mapping;

    private:
        Index_ block_end = 0;

        void fill_mapping() {
            for (size_t i = 0, end = indices.size(); i < end; ++i) {
                mapping[indices[i]] = i;
            }
        }

    public:
        /**
         * @cond
         */
        void set(Index_ i) {
            selection = SubsetSelection::BLOCK;
            block_start = i;
            block_end = i + 1;
            indices.clear();
            mapping.clear();
        }

        void add(Index_ i) {
            if (selection == SubsetSelection::FULL) {
                return;
            }

            if (selection == SubsetSelection::BLOCK) {
                if (i == block_end) {
                    block_end = i + 1;
                    return;

                } else if (i + 1 == block_start) {
                    block_start = i;
                    return;

                } else if (i >= block_start && i < block_end) {
                    return;
                }

                selection = SubsetSelection::INDEX;
                indices.resize(block_end - block_start);
                std::iota(indices.begin(), indices.end(), block_start);
                fill_mapping();
            }

            if (mapping.find(i) == mapping.end()) {
                mapping[i] = indices.size();
                indices.push_back(i);
            }
        }

        void finalize() {
            if (selection == SubsetSelection::BLOCK) {
                block_length = block_end - block_start;
            } else if (selection == SubsetSelection::INDEX) {
                if (!std::is_sorted(indices.begin(), indices.end())) {
                    std::sort(indices.begin(), indices.end());
                    fill_mapping();
                }
            }
        }

        void swap(SubsetDetails& x) {
            std::swap(selection, x.selection);
            std::swap(block_start, x.block_start);
            std::swap(block_end, x.block_end);
            std::swap(block_length, x.block_length);
            indices.swap(x.indices);
            mapping.swap(x.mapping);
        }
        /**
         * @endcond
         */
    };

    /**
     * @brief A cached slab.
     */
    struct CachedSlab {
        /**
         * Contents of the slab.
         */
        Slab_ contents;

        /**
         * Subset of rows/columns to extract from this slab.
         */
        SubsetDetails subset;

    public:
        /**
         * @cond
         */
        CachedSlab(Slab_ c) : contents(std::move(c)) {}
        /**
         * @endcond
         */
    };

private:
    std::list<CachedSlab> slab_cache, tmp_cache, free_cache;

    typedef typename std::list<CachedSlab>::iterator cache_iterator;
    std::unordered_map<Id_, std::pair<Index_, cache_iterator> > slab_exists, past_exists;

    std::vector<std::pair<Index_, Index_> > predictions_made, next_predictions_made;
    size_t predictions_fulfilled = 0;

    std::vector<CachedSlab*> slab_pointers, next_slab_pointers;
    std::vector<std::pair<Index_, cache_iterator*> > unassigned_slabs;
    std::vector<std::pair<Id_, Index_> > slabs_to_populate, next_slabs_to_populate; 

    std::vector<SubsetDetails> next_subset;

public:
    /**
     * @param oracle Pointer to an `tatami::Oracle` to be used for predictions.
     * @param per_iteration Maximum number of predictions to make per iteration.
     * @param num_slabs Maximum number of slabs to store.
     */
    SubsettedOracleSlabCache(std::unique_ptr<tatami::Oracle<Index_> > oracle, size_t per_iteration, size_t num_slabs) :
        prediction_stream(std::move(oracle)), 
        max_predictions(per_iteration),
        max_slabs(num_slabs)
    {
        slab_exists.reserve(max_slabs);
        past_exists.reserve(max_slabs);

        predictions_made.reserve(max_predictions);
        next_predictions_made.reserve(max_predictions);

        slab_pointers.reserve(max_slabs);
        next_slab_pointers.reserve(max_slabs);

        unassigned_slabs.reserve(max_slabs);

        slabs_to_populate.reserve(max_slabs);
        next_slabs_to_populate.reserve(max_slabs);

        next_subset.resize(max_slabs);
    } 

    /**
     * @cond
     */
    // For testing only.
    SubsettedOracleSlabCache() = default;
    /**
     * @endcond
     */

private:
    std::pair<const CachedSlab*, Index_> fetch(size_t i) const {
        const auto& current = predictions_made[i];
        return std::pair<const CachedSlab*, Index_>(slab_pointers[current.first], current.second);
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
     * (2) `slab_data` is a `std::vector<CachedSlab*>&` containing pointers to the cached slab contents to be populated.
     * This function should iterate over the `slabs_in_need` and populate the corresponding entries in `slab_data`,
     * possibly using information in `CachedSlab::subset` to extract only the desired subset of each slab.
     *
     * @return Pair containing (1) a pointer to a cached slab and (2) the index of the next predicted row/column inside the retrieved slab.
     */
    template<class Ifunction_, class Cfunction_, class Pfunction_>
    std::pair<const CachedSlab*, Index_> next(Ifunction_ identify, Cfunction_ create, Pfunction_ populate) {
        if (predictions_made.size() > predictions_fulfilled) {
            return fetch(predictions_fulfilled++);
        }

        if (!next_predictions_made.empty()) {
            predictions_made.swap(next_predictions_made);
            next_predictions_made.clear();

            slab_pointers.swap(next_slab_pointers);
            next_slab_pointers.clear();

            slabs_to_populate.swap(next_slabs_to_populate);
            next_slabs_to_populate.clear();

            // Creating slabs needed for the current prediction round,
            // based on the predictions made in the last round.
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
                *(last.second) = it; // Remember this is a pointer to an iterator, so this assignment changes the iterator in the map.

                unassigned_slabs.pop_back();
            }

            while (!tmp_cache.empty()) {
                free_cache.splice(free_cache.end(), tmp_cache, tmp_cache.begin());
            }

            // Updating subsets for all to-be-populated slabs.
            for (const auto& x : slabs_to_populate) {
                next_subset[x.second].swap(slab_pointers[x.second]->subset);
            }

        } else {
            // This is the first run, so we can freely allocate here.
            size_t used = 0;
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
                    (it->second).second->subset.add(curindex);

                } else if (used < max_slabs) {
                    slab_cache.push_back(CachedSlab(create()));
                    auto sIt = slab_cache.end();
                    --sIt;

                    slab_exists[curslab] = std::make_pair(used, sIt);
                    slabs_to_populate.emplace_back(curslab, used);
                    slab_pointers.push_back(&(*sIt));
                    sIt->subset.set(curindex);

                    predictions_made.emplace_back(used, curindex);
                    ++used;

                } else {
                    prediction_stream.back();
                    break;
                }
            }
        }

        // Now filling up the next round of predictions. Note that this will 
        // change the std::list in which each iterator belongs, along with
        // the various *_exists maps. This is okay as the list iterators and
        // maps are no longer needed for the _current_ prediction round.
        // Only 'slab_pointers' and 'slabs_to_populate' are needed, and these
        // are untouched by the fiddling for the next iteration round.
        size_t used = 0;

        tmp_cache.swap(slab_cache);
        past_exists.swap(slab_exists);
        slab_exists.clear();

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
                auto pos = (it->second).first;
                next_predictions_made.emplace_back(pos, curindex);
                next_subset[pos].add(curindex);

            } else if (used < max_slabs) {
                auto past = past_exists.find(curslab);
                if (past != past_exists.end()) {
                    auto sIt = (past->second).second;
                    slab_cache.splice(slab_cache.end(), tmp_cache, sIt);
                    next_slab_pointers.push_back(&(*sIt));
                    slab_exists[curslab] = std::make_pair(used, sIt);

                    // If we detect that a slab is used in the current and next
                    // prediction rounds, we need to set its subset to FULL
                    // because we don't know whether future iteration rounds
                    // might need even more indices from this slab. The only
                    // way to ensure that this slab is re-usable in subsequent
                    // rounds is to extract it in its entirety.
                    sIt->subset.selection = SubsetSelection::FULL;

                } else {
                    if (free_cache.empty()) {
                        // We might be able to recycle an existing slab from tmp_cache 
                        // to populate 'curslab'... but we don't know if we can do so at
                        // this moment, as those slabs might be needed by later predictions.
                        // So we just defer the creation of a new slab until we've run 
                        // through the set of predictions for this round.
                        auto ins = slab_exists.insert(std::make_pair(curslab, std::make_pair(used, slab_cache.end())));
                        unassigned_slabs.emplace_back(used, &(ins.first->second.second));
                        next_slab_pointers.push_back(NULL);
                        next_subset[used].set(curindex);

                    } else {
                        auto sIt = free_cache.begin();
                        slab_cache.splice(slab_cache.end(), free_cache, sIt);
                        next_slab_pointers.push_back(&(*sIt));
                        slab_exists[curslab] = std::make_pair(used, sIt);
                        next_subset[used].set(curindex);
                    }

                    next_slabs_to_populate.emplace_back(curslab, used);
                }

                next_predictions_made.emplace_back(used, curindex);
                ++used;

            } else {
                prediction_stream.back();
                break;
            }
        }

        // Only populating after the next round of predictions. This is necessary
        // to ensure that subsequent prediction rounds don't need the current slabs;
        // if they don't, we can just extract the specified subset, otherwise we 
        // have to extract the FULL slab to account for potential future use.
        if (!slabs_to_populate.empty()) {
            for (auto& x : slabs_to_populate) {
                slab_pointers[x.second]->subset.finalize();
            }
            populate(slabs_to_populate, slab_pointers);
        }

        // Well, because we just used one.
        predictions_fulfilled = 1;
        return fetch(0);
    }
};

}

#endif
