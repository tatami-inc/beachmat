#ifndef TATAMI_CHUNKED_ORACULAR_VARIABLE_SLAB_CACHE_HPP
#define TATAMI_CHUNKED_ORACULAR_VARIABLE_SLAB_CACHE_HPP

#include <unordered_map>
#include <vector>
#include <list>
#include <type_traits>
#include "tatami/tatami.hpp"

/**
 * @file OracularVariableSlabCache.hpp
 * @brief Oracle-aware cache for variable-size slabs.
 */

namespace tatami_chunked {

/**
 * @brief Oracle-aware cache for variable-size slabs.
 *
 * @tparam Id_ Type of slab identifier, typically integer.
 * @tparam Index_ Type of row/column index produced by the oracle.
 * @tparam Size_ Numeric type for the slab size.
 * @tparam Slab_ Class for a single slab.
 *
 * Implement an oracle-aware cache for variable-size slabs.
 * Each slab is defined as the set of chunks required to read an element of the target dimension (or a contiguous block/indexed subset thereof) from a `tatami::Matrix`.
 * This cache is similar to `OracularSlabCache` but enables improved cache utilization when the slabs vary in size.
 * For example, the number of non-zero entries in a sparse matrix might vary between slabs,
 * so the cache could be optimized to fit more slabs into memory when they have fewer non-zeros.
 *
 * The size of each slab is defined by `Size_`, which can be any non-negative measure of slab size.
 * This could be the number of non-zero elements, or the number of dimension elements, or the size of the slab in bytes, etc.,
 * as long as its interpretation is consistent between slabs and with the `max_size` used in the constructor.
 * Users can also differentiate between the estimated and actual size of the slab, if the latter is not known until after the slab has been loaded into memory,
 * e.g., the number of non-zero entries in a file-backed sparse matrix.
 *
 * When implementing `Slab_`,  we generally suggest using a common memory pool that is referenced by each `Slab_` instance.
 * This guarantees that the actual cache size does not exceed the limit associated with `max_size` when `Slab_` instances are re-used for different slabs.
 * (Otherwise, if each `Slab_` allocates its own memory, re-use of an instance may cause its allocation to increase to the size of the largest enmy_countered slab.)
 * Callers may need to occasionally defragment the pool to ensure that enough memory is available for loading new slabs.
 */
template<typename Id_, typename Index_, class Slab_, typename Size_> 
class OracularVariableSlabCache {
private:
    std::shared_ptr<const tatami::Oracle<Index_> > my_oracle;
    size_t my_total;
    size_t my_counter = 0;

    Index_ my_last_slab_id = 0;
    size_t my_last_slab_num = -1;

    Size_ my_max_size, my_used_size = 0;
    std::vector<Slab_> my_all_slabs;

    // We need to hold an offset into 'my_all_slabs' rather than a pointer, as
    // 'my_all_slabs' might be reallocated upon addition of new slabs, given that
    // we don't know the maximum number of slabs ahead of time.
    std::unordered_map<Id_, size_t> my_current_cache, my_future_cache;
    std::vector<std::pair<Id_, size_t> > my_to_populate, my_to_reuse;
    std::vector<Id_> my_in_need;
    std::vector<size_t> my_free_pool;
    size_t my_refresh_point = 0;

public:
    /**
     * @param oracle Pointer to an `tatami::Oracle` to be used for predictions.
     * @param max_size Total size of all slabs to store in the cache.
     * This may be zero, in which case no caching should be performed.
     */
    OracularVariableSlabCache(std::shared_ptr<const tatami::Oracle<Index_> > oracle, size_t max_size) : 
        my_oracle(std::move(oracle)), 
        my_total(my_oracle->total()),
        my_max_size(max_size) 
    {} 

    /**
     * Deleted as the cache holds persistent pointers.
     */
    OracularVariableSlabCache(const OracularVariableSlabCache&) = delete;

    /**
     * Deleted as the cache holds persistent pointers.
     */
    OracularVariableSlabCache& operator=(const OracularVariableSlabCache&) = delete;

    /**
     * @cond
     */
    // Move operators are still okay as pointers still point to the moved vectors.
    // see https://stackoverflow.com/questions/43988553/stdvector-stdmove-and-pointer-invalidation.
    OracularVariableSlabCache& operator=(OracularVariableSlabCache&&) = default;
    OracularVariableSlabCache(OracularVariableSlabCache&&) = default;

    // Might as well define this.
    ~OracularVariableSlabCache() = default;
    /**
     * @endcond
     */

public:
    /**
     * This method is intended to be called when `max_size = 0`, to provide callers with the oracle predictions for non-cached extraction of data.
     * Calls to this method should not be intermingled with calls to its overload below; the latter should only be called when `max_size > 0`.
     *
     * @return The next prediction from the oracle.
     */
    Index_ next() {
        return my_oracle->get(my_counter++);
    }

public:
    /**
     * Fetch the next slab according to the stream of predictions provided by the `tatami::Oracle`.
     * This method should only be called if `max_size > 0` in the constructor; otherwise, no slabs are actually available and cannot be returned.
     *
     * @tparam Ifunction_ Function to identify the slab containing each predicted row/column.
     * @tparam Efunction_ Function to compute the estimated size of a slab.
     * @tparam Afunction_ Function to compute the actual size of a slab.
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
     * @param estimated_size Function that accepts `j`, an `Id_` containing the slab identifier.
     * It should return the size of the slab as a non-negative `Size_`.
     * @param actual_size Function that accepts `j`, an `Id_` containing the slab identifier; and `slab`, a populated `const Slab_&` instance corresponding to `j`.
     * It should return the actual size of the slab as a non-negative `Size_` that is no greater than `estimated_size(j)`.
     * @param create Function that accepts no arguments and returns a `Slab_` object with sufficient memory to hold a slab's contents when used in `populate()`.
     * This may also return a default-constructed `Slab_` object if the allocation is done dynamically per slab in `populate()`.
     * @param populate Function that accepts three arguments - `my_to_populate`, `my_to_reuse` and `my_all_slabs`.
     * - The `my_to_populate` argument is a `std::vector<std::pair<Id_, size_t> >&` specifying the slabs to be populated.
     *   The first `Id_` element of each pair contains the slab identifier, i.e., the first element returned by the `identify` function.
     *   The second `size_t` element specifies the entry of `my_all_slabs` containing the corresponding `Slab_` instance, as returned by `create()`.
     *   This argument can be modified in any manner.
     *   It is not guaranteed to be sorted.
     * - The `my_to_reuse` argument is a `std::vector<std::pair<Id_, size_t> >&` specifying the cached slabs that were re-used in the upcoming set of predictions.
     *   The elements of each pair are interpreted in the same manner as `my_to_populate`. 
     *   This argument can be modified in any manner.
     *   It is not guaranteed to be sorted.
     * - The `my_all_slabs` argument is a `std::vector<Slab_>&` containing all slabs in the cache.
     *   This may include instances that are not referenced by `my_to_populate` or `my_to_reuse`.
     *   Each element of this argument can be modified but the length should not change.
     * .
     * The `populate` function should iterate over `my_to_populate` and fill each `Slab_` with the contents of the corresponding slab.
     * Optionally, callers may use `my_to_reuse` to defragment the already-in-use parts of the cache, in order to free up enough space for new data from `my_to_populate`.
     *
     * @return Pair containing (1) a pointer to a slab's contents and (2) the index of the next predicted row/column inside the retrieved slab.
     */
    template<class Ifunction_, class Efunction_, class Afunction_, class Cfunction_, class Pfunction_>
    std::pair<const Slab_*, Index_> next(Ifunction_ identify, Efunction_ estimated_size, Afunction_ actual_size, Cfunction_ create, Pfunction_ populate) {
        Index_ index = this->next(); 
        auto slab_info = identify(index);
        if (slab_info.first == my_last_slab_id && my_last_slab_num != static_cast<size_t>(-1)) {
            return std::make_pair(my_all_slabs.data() + my_last_slab_num, slab_info.second);
        }
        my_last_slab_id = slab_info.first;

        // Updating the cache if we hit the refresh point.
        if (my_counter - 1 == my_refresh_point) {
            // Note that, for any given populate cycle, the first prediction's
            // slab cannot already be in the cache, otherwise it would have
            // incorporated into the previous cycle. So we can skip some code.
            my_used_size = estimated_size(slab_info.first);
            requisition_new_slab(slab_info.first);

            auto last_future_slab_id = slab_info.first;
            while (++my_refresh_point < my_total) {
                auto future_index = my_oracle->get(my_refresh_point);
                auto future_slab_info = identify(future_index);
                if (last_future_slab_id == future_slab_info.first) {
                    continue;
                }

                last_future_slab_id = future_slab_info.first;
                if (my_future_cache.find(future_slab_info.first) != my_future_cache.end()) {
                    continue;
                }

                auto ccIt = my_current_cache.find(future_slab_info.first);
                if (ccIt != my_current_cache.end()) {
                    size_t slab_num = ccIt->second;
                    auto candidate = my_used_size + actual_size(future_slab_info.first, my_all_slabs[slab_num]);
                    if (candidate > my_max_size) {
                        break;
                    } 
                    my_used_size = candidate;
                    my_future_cache[future_slab_info.first] = slab_num;
                    my_to_reuse.emplace_back(future_slab_info.first, slab_num);
                    my_current_cache.erase(ccIt);
                } else {
                    auto candidate = my_used_size + estimated_size(future_slab_info.first);
                    if (candidate > my_max_size) {
                        break;
                    } 
                    my_used_size = candidate;
                    requisition_new_slab(future_slab_info.first);
                }
            }

            auto cIt = my_current_cache.begin();
            for (auto a : my_in_need) {
                if (cIt != my_current_cache.end()) {
                    size_t slab_num = cIt->second;
                    my_to_populate.emplace_back(a, slab_num);
                    my_future_cache[a] = slab_num;
                    ++cIt;
                } else {
                    size_t slab_num = my_all_slabs.size();
                    my_all_slabs.push_back(create());
                    my_to_populate.emplace_back(a, slab_num);
                    my_future_cache[a] = slab_num;
                }
            }
            my_in_need.clear();

            for (; cIt != my_current_cache.end(); ++cIt) {
                my_free_pool.emplace_back(cIt->second);
            }

            populate(my_to_populate, my_to_reuse, my_all_slabs);
            my_to_populate.clear();
            my_to_reuse.clear();

            my_current_cache.clear();
            my_current_cache.swap(my_future_cache);
        }

        // We know it must exist, so no need to check ccIt's validity.
        auto ccIt = my_current_cache.find(slab_info.first);
        my_last_slab_num = ccIt->second;
        return std::make_pair(my_all_slabs.data() + my_last_slab_num, slab_info.second);
    }

private:
    void requisition_new_slab(Id_ slab_id) {
        if (!my_free_pool.empty()) {
            auto slab_num = my_free_pool.back();
            my_future_cache[slab_id] = slab_num;
            my_free_pool.pop_back();
            my_to_populate.emplace_back(slab_id, slab_num);
        } else {
            my_future_cache[slab_id] = 0;
            my_in_need.push_back(slab_id);
        }
    }

public:
    /**
     * @return Maximum total size of the cache.
     * This is the same as the `max_size` used in the constructor.
     */
    size_t get_max_size() const {
        return my_max_size;
    }

    /**
     * @return Current usage across all slabs in the cache.
     * This should be interpreted as an upper bound on usage if there is a difference between estimated and actual slab sizes.
     */
    size_t get_used_size() const {
        return my_used_size;
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
