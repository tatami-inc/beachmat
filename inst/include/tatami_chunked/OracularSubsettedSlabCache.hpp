#ifndef TATAMI_CHUNKED_SUBSETTED_ORACLE_SLAB_CACHE_HPP
#define TATAMI_CHUNKED_SUBSETTED_ORACLE_SLAB_CACHE_HPP

#include <unordered_map>
#include <vector>
#include <list>
#include "tatami/tatami.hpp"

/**
 * @file OracularSubsettedSlabCache.hpp
 * @brief Create a oracle-aware cache with subsets.
 */

namespace tatami_chunked {

/**
 * Type of selection on the target dimension.
 * Used to determine the subsets to extract in `OracularSubsettedSlabCache`.
 * 
 * - `FULL`: all rows/columns. 
 * - `BLOCK`: a contiguous block of rows/columns. 
 * - `INDEX`: an indexed subset of rows/columns.
 */
enum class OracularSubsettedSlabCacheSelectionType : char { FULL, BLOCK, INDEX };

/**
 * @brief Details on the subset to extract in `OracularSubsettedSlabCache`.
 * @tparam Index_ Type of row/column index produced by the oracle.
 */
template<typename Index_>
struct OracularSubsettedSlabCacheSelectionDetails {
    /**
     * Type of subset to extract from the target dimension of the slab.
     */
    OracularSubsettedSlabCacheSelectionType selection;

    /**
     * Row/column index representing the start of the contiguous block of the target dimension to be extracted from the slab.
     * Only used if `selection` is set to `OracularSubsettedSlabCacheSelectionType::BLOCK`.
     *
     * Note that the index is relative to the start of the slab, not to the matrix containing the slab,
     * i.e., if the slab consists of rows 10-20 and we want to extract row 11, this will be reported here as an index of 1.
     */
    Index_ block_start;

    /**
     * Length of the contiguous block of the target dimension to be extracted from the slab.
     * Only used if `selection` is set to `OracularSubsettedSlabCacheSelectionType::BLOCK`.
     */
    Index_ block_length;

    /**
     * Row/column index representing one-past-the-end of the contiguous block of the target dimension to be extracted from the slab.
     * Only used if `selection` is set to `OracularSubsettedSlabCacheSelectionType::BLOCK`.
     * This is also equal to `block_start + block_length`.
     */
    Index_ block_end;

    /**
     * Indices of the target dimension to be extracted from the slab.
     * Guaranteed to be sorted and unique.
     * Only used if `selection` is set to `OracularSubsettedSlabCacheSelectionType::INDEX`.
     *
     * Note that all indices is relative to the start of the slab, not to the matrix containing the slab,
     * i.e., if the slab consists of rows 10-20 and we want to extract row 11, this will be reported here as an index of 1.
     */
    std::vector<Index_> indices;

    /**
     * Mapping of indices-to-be-extracted to their positions inside `indices`.
     * All values of `indices` are present as keys here where `mapping[indices[i]] = i`.
     * Only used if `selection` is set to `OracularSubsettedSlabCacheSelectionType::INDEX`.
     */
    std::unordered_map<Index_, Index_> mapping;
};

/**
 * @cond
 */
namespace OracularSubsettedSlabCache_internals {

// We put these functions in here as struct{} usage policy forbids methods
// (outside of the constructor). The Details class should be a passive data
// carrier only.

template<typename Index_>
void fill_mapping_in_details(OracularSubsettedSlabCacheSelectionDetails<Index_>& details) {
    for (size_t i = 0, end = details.indices.size(); i < end; ++i) {
        details.mapping[details.indices[i]] = i;
    }
}

template<typename Index_>
void set_details(OracularSubsettedSlabCacheSelectionDetails<Index_>& details, Index_ i) {
    details.selection = OracularSubsettedSlabCacheSelectionType::BLOCK;
    details.block_start = i;
    details.block_end = i + 1;
    details.indices.clear();
    details.mapping.clear();
}

template<typename Index_>
void add_to_details(OracularSubsettedSlabCacheSelectionDetails<Index_>& details, Index_ i) {
    if (details.selection == OracularSubsettedSlabCacheSelectionType::FULL) {
        return;
    }

    if (details.selection == OracularSubsettedSlabCacheSelectionType::BLOCK) {
        if (i == details.block_end) {
            details.block_end = i + 1;
            return;

        } else if (i + 1 == details.block_start) {
            details.block_start = i;
            return;

        } else if (i >= details.block_start && i < details.block_end) {
            return;
        }

        details.selection = OracularSubsettedSlabCacheSelectionType::INDEX;
        details.indices.resize(details.block_end - details.block_start);
        std::iota(details.indices.begin(), details.indices.end(), details.block_start);
        fill_mapping_in_details(details);
    }

    if (details.mapping.find(i) == details.mapping.end()) {
        details.mapping[i] = details.indices.size();
        details.indices.push_back(i);
    }
}

template<typename Index_>
void finalize_details(OracularSubsettedSlabCacheSelectionDetails<Index_>& details) {
    if (details.selection == OracularSubsettedSlabCacheSelectionType::BLOCK) {
        details.block_length = details.block_end - details.block_start;
    } else if (details.selection == OracularSubsettedSlabCacheSelectionType::INDEX) {
        if (!std::is_sorted(details.indices.begin(), details.indices.end())) {
            std::sort(details.indices.begin(), details.indices.end());
            fill_mapping_in_details(details);
        }
    }
}

}
/**
 * @endcond
 */


/**
 * @brief Oracle-aware cache for slabs, plus subsets.
 *
 * @tparam Id_ Type of slab identifier, typically integer.
 * @tparam Index_ Type of row/column index produced by the oracle.
 * @tparam Slab_ Class for a single slab.
 *
 * Implement an oracle-aware cache for slab subsets.
 * Each slab is defined as the set of chunks required to read an element of the target dimension (or a contiguous block/indexed subset thereof) from a `tatami::Matrix`.
 * This cache is similar to the `OracularSlabCache` except that it remembers the subset of elements on the target dimension that were requested for each slab.
 * Slab extractors can use this information to optimize slab loading by ignoring unneeded elements. 
 */
template<typename Id_, typename Index_, class Slab_> 
class OracularSubsettedSlabCache {
private:
    std::shared_ptr<const tatami::Oracle<Index_> > my_oracle;
    size_t my_total;
    size_t my_counter = 0;

    Index_ my_last_slab_id = 0;
    Slab_* my_last_slab = NULL;

    size_t my_max_slabs;
    std::vector<Slab_> my_all_slabs;
    std::unordered_map<Id_, Slab_*> my_current_cache, my_future_cache;

    std::vector<OracularSubsettedSlabCacheSelectionDetails<Index_> > my_all_subset_details;
    std::vector<OracularSubsettedSlabCacheSelectionDetails<Index_>*> my_free_subset_details;
    std::unordered_map<Id_, OracularSubsettedSlabCacheSelectionDetails<Index_>*> my_close_future_subset_cache, my_far_future_subset_cache;
    size_t my_close_refresh_point = 0;
    size_t my_far_refresh_point = 0;
    Id_ my_far_slab_id;
    Index_ my_far_slab_offset;

    std::vector<std::pair<Id_, OracularSubsettedSlabCacheSelectionDetails<Index_>*> > my_to_reassign;
    std::vector<std::tuple<Id_, Slab_*, const OracularSubsettedSlabCacheSelectionDetails<Index_>*> > my_to_populate;

public:
    /**
     * @param oracle Pointer to an `tatami::Oracle` to be used for predictions.
     * @param max_slabs Maximum number of slabs to store.
     */
    OracularSubsettedSlabCache(std::shared_ptr<const tatami::Oracle<Index_> > oracle, size_t max_slabs) :
        my_oracle(std::move(oracle)), 
        my_total(my_oracle->total()),
        my_max_slabs(max_slabs)
    {
        my_all_slabs.reserve(max_slabs);
        my_current_cache.reserve(max_slabs);
        my_future_cache.reserve(max_slabs);
        my_close_future_subset_cache.reserve(max_slabs);
        my_far_future_subset_cache.reserve(max_slabs);

        my_all_subset_details.resize(max_slabs * 2);
        for (auto& as : my_all_subset_details) {
            my_free_subset_details.push_back(&as);
        }
    }

    /**
     * Deleted as the cache holds persistent pointers.
     */
    OracularSubsettedSlabCache(const OracularSubsettedSlabCache&) = delete;

    /**
     * Deleted as the cache holds persistent pointers.
     */
    OracularSubsettedSlabCache& operator=(const OracularSubsettedSlabCache&) = delete;

    /**
     * @cond
     */
    // Move operators are still okay as pointers still point to the moved vectors,
    // see https://stackoverflow.com/questions/43988553/stdvector-stdmove-and-pointer-invalidation.
    OracularSubsettedSlabCache(OracularSubsettedSlabCache&&) = delete;
    OracularSubsettedSlabCache& operator=(OracularSubsettedSlabCache&&) = delete;

    // Might as well define this.
    ~OracularSubsettedSlabCache() = default;
    /**
     * @endcond
     */

public:
    /**
     * This method is intended to be called when `num_slabs = 0`, to provide callers with the oracle predictions for non-cached extraction of data.
     * Calls to this method should not be intermingled with calls to its overload below; the latter should only be called when `max_slabs > 0`.
     *
     * @return The next prediction from the oracle.
     */
    Index_ next() {
        return my_oracle->get(my_counter++);
    }

public:
    /**
     * Fetch the next slab according to the stream of predictions provided by the `tatami::Oracle`.
     * This method should only be called if `num_slabs > 0` in the constructor; otherwise, no slabs are actually available and cannot be returned.
     *
     * @tparam Ifunction_ Function to identify the slab containing each predicted row/column.
     * @tparam Cfunction_ Function to create a new slab.
     * @tparam Pfunction_ Function to populate zero, one or more slabs with their contents.
     *
     * @param identify Function that accepts an `i`, an `Index_` containing the predicted index of a single element on the target dimension.
     * This should return a pair containing:
     * 1. An `Id_`, the identifier of the slab containing `i`.
     *    This is typically defined as the index of the slab on the target dimension.
     *    For example, if each chunk takes up 10 rows, attempting to access row 21 would require retrieval of slab 2.
     * 2. An `Index_`, the index of row/column `i` inside that slab.
     *    For example, if each chunk takes up 10 rows, attempting to access row 21 would yield an offset of 1.
     * @param create Function that accepts no arguments and returns a `Slab_` object with sufficient memory to hold a slab's contents when used in `populate()`.
     * This may also return a default-constructed `Slab_` object if the allocation is done dynamically per slab in `populate()`.
     * @param populate Function that accepts a `std::vector<std::tuple<Id_, Slab_*, const OracularSubsettedSlabCacheSelectionDetails<Index_>*> >&` specifying the slabs to be populated.
     * The first `Id_` element of each tuple contains the slab identifier, i.e., the first element returned by the `identify` function.
     * The second `Slab_*` element specifies the object which to store the contents of the slab.
     * The third `OracularSubsettedSlabCacheSelectionDetails<Index_>*` element contains information about the desired subset of elements on the target dimension of the slab.
     * This function should iterate over the vector and populate the desired subset of each slab.
     * Note that the vector is not guaranteed to be sorted. 
     *
     * @return Pair containing (1) a pointer to a cached slab and (2) the index of the next predicted row/column inside the retrieved slab.
     */
    template<class Ifunction_, class Cfunction_, class Pfunction_>
    std::pair<const Slab_*, Index_> next(Ifunction_ identify, Cfunction_ create, Pfunction_ populate) {
        Index_ index = this->next(); 
        auto slab_info = identify(index);
        if (slab_info.first == my_last_slab_id && my_last_slab) {
            return std::make_pair(my_last_slab, slab_info.second);
        }
        my_last_slab_id = slab_info.first;

        // Updating the cache if we hit the refresh point.
        if (my_counter - 1 == my_close_refresh_point) {
            if (my_all_slabs.empty()) {
                // This section only runs once, at the start, to populate the my_close_future_subset_cache.
                requisition_subset_close(slab_info.first, slab_info.second);
                size_t used_slabs = 1;

                while (++my_close_refresh_point < my_total) {
                    auto future_index = my_oracle->get(my_close_refresh_point);
                    auto future_slab_info = identify(future_index);
                    auto cfcIt = my_close_future_subset_cache.find(future_slab_info.first);
                    if (cfcIt != my_close_future_subset_cache.end()) {
                        OracularSubsettedSlabCache_internals::add_to_details(*(cfcIt->second), future_slab_info.second);
                    } else if (used_slabs < my_max_slabs) {
                        requisition_subset_close(future_slab_info.first, future_slab_info.second);
                        ++used_slabs;
                    } else {
                        my_far_slab_id = future_slab_info.first;
                        my_far_slab_offset = future_slab_info.second;
                        break;
                    }
                }

                my_far_refresh_point = my_close_refresh_point;
            } else {
                my_close_refresh_point = my_far_refresh_point;
            }

            // Populating the far future cache. 
            if (my_far_refresh_point < my_total) {
                requisition_subset_far(my_far_slab_id, my_far_slab_offset);
                size_t used_slabs = 1;

                while (++my_far_refresh_point < my_total) {
                    auto future_index = my_oracle->get(my_far_refresh_point);
                    auto future_slab_info = identify(future_index);
                    auto ffcIt = my_far_future_subset_cache.find(future_slab_info.first);
                    if (ffcIt != my_far_future_subset_cache.end()) {
                        OracularSubsettedSlabCache_internals::add_to_details(*(ffcIt->second), future_slab_info.second);
                    } else if (used_slabs < my_max_slabs) {
                        requisition_subset_far(future_slab_info.first, future_slab_info.second);
                        ++used_slabs;
                    } else {
                        my_far_slab_id = future_slab_info.first;
                        my_far_slab_offset = future_slab_info.second;
                        break;
                    }
                }
            }

            // Reusing slabs from my_current_cache; these should all have FULL selections already.
            for (auto& cf : my_close_future_subset_cache) {
                auto cIt = my_current_cache.find(cf.first);
                if (cIt == my_current_cache.end()) {
                    my_to_reassign.emplace_back(cf.first, cf.second);
                } else {
                    my_future_cache[cf.first] = cIt->second;
                    my_current_cache.erase(cIt);
                }
            }

            // Creating new slabs for everything that's left.
            auto cIt = my_current_cache.begin();
            for (auto a : my_to_reassign) {
                Slab_* slab_ptr;
                if (cIt == my_current_cache.end()) {
                    my_all_slabs.emplace_back(create());
                    slab_ptr = &(my_all_slabs.back());
                } else {
                    slab_ptr = cIt->second;
                    ++cIt;
                }
                my_future_cache[a.first] = slab_ptr;
                OracularSubsettedSlabCache_internals::finalize_details(*(a.second));
                my_to_populate.emplace_back(a.first, slab_ptr, a.second);
            }
            my_to_reassign.clear();

            populate(my_to_populate);
            my_to_populate.clear();

            // We always fill my_future_cache to the brim so every entry of
            // my_all_slabs should be referenced by a pointer in
            // my_future_cache.  There shouldn't be any free cache entries
            // remaining in my_current_cache i.e., at this point, cIt should
            // equal my_current_cache.end(), as we transferred everything to
            // my_future_cache. Thus it is safe to clear my_current_cache
            // without worrying about leaking memory. The only exception is if
            // we run out of predictions, in which case it doesn't matter.
            my_current_cache.clear();
            my_current_cache.swap(my_future_cache);

            // Putting the no-longer-used subset pointers back in the free pool
            // before we swap the close and far futures.
            for (auto& cfc : my_close_future_subset_cache) {
                my_free_subset_details.push_back(cfc.second);
            }
            my_close_future_subset_cache.clear();
            my_close_future_subset_cache.swap(my_far_future_subset_cache);
        }

        // We know it must exist, so no need to check ccIt's validity.
        auto ccIt = my_current_cache.find(slab_info.first);
        my_last_slab = ccIt->second;
        return std::make_pair(my_last_slab, slab_info.second);
    }

private:
    void requisition_subset_close(Id_ slab_id, Index_ slab_offset) {
        auto selected = my_free_subset_details.back();
        OracularSubsettedSlabCache_internals::set_details(*selected, slab_offset);
        my_close_future_subset_cache[slab_id] = selected;
        my_free_subset_details.pop_back();
    }

    void requisition_subset_far(Id_ slab_id, Index_ slab_offset) {
        auto selected = my_free_subset_details.back();
        OracularSubsettedSlabCache_internals::set_details(*selected, slab_offset);
        my_far_future_subset_cache[slab_id] = selected;
        my_free_subset_details.pop_back();

        // If a slab is still being used in the far future, it might continue
        // to be used in an even further future, in which case we need to do a
        // FULL extraction just to be safe.
        auto cfcIt = my_close_future_subset_cache.find(slab_id);
        if (cfcIt != my_close_future_subset_cache.end()) {
            selected->selection = OracularSubsettedSlabCacheSelectionType::FULL;
            cfcIt->second->selection = OracularSubsettedSlabCacheSelectionType::FULL;
        }
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
        return my_current_cache.size();
    }
};

}

#endif
