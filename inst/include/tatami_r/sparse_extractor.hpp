#ifndef TATAMI_R_SPARSE_EXTRACTOR_HPP
#define TATAMI_R_SPARSE_EXTRACTOR_HPP

#include "Rcpp.h"
#include "tatami/tatami.hpp"
#include "tatami_chunked/tatami_chunked.hpp"
#include "sparse_matrix.hpp"

#include <vector>
#include <stdexcept>

namespace tatami_r {

namespace UnknownMatrix_internal {

/********************
 *** Core classes ***
 ********************/

template<bool oracle_, typename Index_, typename CachedValue_, typename CachedIndex_>
class SoloSparseCore {
public:
    SoloSparseCore(
        const Rcpp::RObject& matrix, 
        const Rcpp::Function& sparse_extractor,
        bool row,
        tatami::MaybeOracle<oracle_, Index_> oracle,
        Rcpp::IntegerVector non_target_extract, 
        [[maybe_unused]] Index_ max_target_chunk_length, // provided here for compatibility with the other Sparse*Core classes.
        [[maybe_unused]] const std::vector<Index_>& ticks,
        [[maybe_unused]] const std::vector<Index_>& map,
        [[maybe_unused]] const tatami_chunked::SlabCacheStats& stats,
        bool needs_value,
        bool needs_index) : 
        my_matrix(matrix),
        my_sparse_extractor(sparse_extractor),
        my_extract_args(2),
        my_row(row),
        my_factory(1, non_target_extract.size(), 1, needs_value, needs_index),
        my_solo(my_factory.create()),
        my_oracle(std::move(oracle))
    {
        my_extract_args[static_cast<int>(row)] = non_target_extract;
    }

private:
    const Rcpp::RObject& my_matrix;
    const Rcpp::Function& my_sparse_extractor;
    Rcpp::List my_extract_args;

    bool my_row;

    tatami_chunked::SparseSlabFactory<CachedValue_, CachedIndex_> my_factory;
    typedef typename decltype(my_factory)::Slab Slab;
    Slab my_solo;

    tatami::MaybeOracle<oracle_, Index_> my_oracle;
    typename std::conditional<oracle_, size_t, bool>::type my_counter = 0;

public:
    std::pair<const Slab*, Index_> fetch_raw(Index_ i) {
        if constexpr(oracle_) {
            i = my_oracle->get(my_counter++);
        }
        my_solo.number[0] = 0;

#ifdef TATAMI_R_PARALLELIZE_UNKNOWN 
        // This involves some Rcpp initializations, so we lock it just in case.
        auto& mexec = executor();
        mexec.run([&]() -> void {
#endif

        my_extract_args[static_cast<int>(!my_row)] = Rcpp::IntegerVector::create(i + 1);
        auto obj = my_sparse_extractor(my_matrix, my_extract_args);
        parse_sparse_matrix(obj, my_row, my_solo.values, my_solo.indices, my_solo.number);

#ifdef TATAMI_R_PARALLELIZE_UNKNOWN 
        });
#endif

        return std::make_pair(&my_solo, static_cast<Index_>(0));
    }
};

template<typename Index_, typename CachedValue_, typename CachedIndex_>
class MyopicSparseCore {
public:
    MyopicSparseCore(
        const Rcpp::RObject& matrix, 
        const Rcpp::Function& sparse_extractor,
        bool row,
        [[maybe_unused]] tatami::MaybeOracle<false, Index_> oracle, // provided here for compatibility with the other Sparse*Core classes.
        Rcpp::IntegerVector non_target_extract, 
        Index_ max_target_chunk_length, 
        const std::vector<Index_>& ticks,
        const std::vector<Index_>& map,
        const tatami_chunked::SlabCacheStats& stats,
        bool needs_value,
        bool needs_index) : 
        my_matrix(matrix),
        my_sparse_extractor(sparse_extractor),
        my_extract_args(2),
        my_row(row),
        my_chunk_ticks(ticks),
        my_chunk_map(map),
        my_factory(max_target_chunk_length, non_target_extract.size(), stats, needs_value, needs_index),
        my_cache(stats.max_slabs_in_cache)
    {
        my_extract_args[static_cast<int>(row)] = non_target_extract;
    }

private:
    const Rcpp::RObject& my_matrix;
    const Rcpp::Function& my_sparse_extractor;
    Rcpp::List my_extract_args;

    bool my_row;

    const std::vector<Index_>& my_chunk_ticks;
    const std::vector<Index_>& my_chunk_map;

    tatami_chunked::SparseSlabFactory<CachedValue_, CachedIndex_> my_factory;
    typedef typename decltype(my_factory)::Slab Slab;
    tatami_chunked::LruSlabCache<Index_, Slab> my_cache;

public:
    std::pair<const Slab*, Index_> fetch_raw(Index_ i) {
        auto chosen = my_chunk_map[i];

        const auto& slab = my_cache.find(
            chosen,
            [&]() -> Slab {
                return my_factory.create();
            },
            [&](Index_ id, Slab& cache) {
                auto chunk_start = my_chunk_ticks[id], chunk_end = my_chunk_ticks[id + 1];
                size_t chunk_len = chunk_end - chunk_start;
                std::fill_n(cache.number, chunk_len, 0);

#ifdef TATAMI_R_PARALLELIZE_UNKNOWN 
                // This involves some Rcpp initializations, so we lock it just in case.
                auto& mexec = executor();
                mexec.run([&]() -> void {
#endif

                Rcpp::IntegerVector primary_extract(chunk_len);
                std::iota(primary_extract.begin(), primary_extract.end(), chunk_start + 1);
                my_extract_args[static_cast<int>(!my_row)] = primary_extract;
                auto obj = my_sparse_extractor(my_matrix, my_extract_args);
                parse_sparse_matrix(obj, my_row, cache.values, cache.indices, cache.number);

#ifdef TATAMI_R_PARALLELIZE_UNKNOWN 
                });
#endif
            }
        );

        Index_ offset = i - my_chunk_ticks[chosen];
        return std::make_pair(&slab, offset);
    }
};

template<typename Index_, typename CachedValue_, typename CachedIndex_>
class OracularSparseCore {
public:
    OracularSparseCore(
        const Rcpp::RObject& matrix, 
        const Rcpp::Function& sparse_extractor,
        bool row,
        tatami::MaybeOracle<true, Index_> oracle,
        Rcpp::IntegerVector non_target_extract, 
        Index_ max_target_chunk_length, 
        const std::vector<Index_>& ticks,
        const std::vector<Index_>& map,
        const tatami_chunked::SlabCacheStats& stats,
        bool needs_value,
        bool needs_index) : 
        my_matrix(matrix),
        my_sparse_extractor(sparse_extractor),
        my_extract_args(2),
        my_row(row),
        my_chunk_ticks(ticks),
        my_chunk_map(map),
        my_factory(max_target_chunk_length, non_target_extract.size(), stats, needs_value, needs_index),
        my_cache(std::move(oracle), stats.max_slabs_in_cache),
        my_needs_value(needs_value),
        my_needs_index(needs_index)
    {
        my_extract_args[static_cast<int>(row)] = non_target_extract;
    }

private:
    const Rcpp::RObject& my_matrix;
    const Rcpp::Function& my_sparse_extractor;
    Rcpp::List my_extract_args;

    bool my_row;

    const std::vector<Index_>& my_chunk_ticks;
    const std::vector<Index_>& my_chunk_map;

    tatami_chunked::SparseSlabFactory<CachedValue_, CachedIndex_> my_factory;
    typedef typename decltype(my_factory)::Slab Slab;
    tatami_chunked::OracularSlabCache<Index_, Index_, Slab> my_cache;

    std::vector<CachedValue_*> my_chunk_value_ptrs;
    std::vector<CachedIndex_*> my_chunk_index_ptrs;
    std::vector<CachedIndex_> my_chunk_numbers;

    bool my_needs_value;
    bool my_needs_index;

public:
    std::pair<const Slab*, Index_> fetch_raw(Index_) {
        return my_cache.next(
            [&](Index_ i) -> std::pair<Index_, Index_> {
                auto chosen = my_chunk_map[i];
                return std::make_pair(chosen, static_cast<Index_>(i - my_chunk_ticks[chosen]));
            },
            [&]() -> Slab {
                return my_factory.create();
            },
            [&](std::vector<std::pair<Index_, Slab*> >& to_populate) {
                // Sorting them so that the indices are in order.
                if (!std::is_sorted(to_populate.begin(), to_populate.end(), [&](const std::pair<Index_, Slab*>& left, const std::pair<Index_, Slab*> right) { return left.first < right.first; })) {
                    std::sort(to_populate.begin(), to_populate.end(), [&](const std::pair<Index_, Slab*>& left, const std::pair<Index_, Slab*> right) { return left.first < right.first; });
                }

                if (my_needs_value) {
                    my_chunk_value_ptrs.clear();
                }
                if (my_needs_index) {
                    my_chunk_index_ptrs.clear();
                }

                Index_ total_len = 0;
                for (const auto& p : to_populate) {
                    Index_ chunk_len = my_chunk_ticks[p.first + 1] - my_chunk_ticks[p.first];
                    total_len += chunk_len;
                    if (my_needs_value) {
                        auto vIt = p.second->values.begin();
                        my_chunk_value_ptrs.insert(my_chunk_value_ptrs.end(), vIt, vIt + chunk_len);
                    }
                    if (my_needs_index) {
                        auto iIt = p.second->indices.begin();
                        my_chunk_index_ptrs.insert(my_chunk_index_ptrs.end(), iIt, iIt + chunk_len);
                    }
                }

                my_chunk_numbers.clear();
                my_chunk_numbers.resize(total_len);

#ifdef TATAMI_R_PARALLELIZE_UNKNOWN 
                // This involves some Rcpp initializations, so we lock it just in case.
                auto& mexec = executor();
                mexec.run([&]() -> void {
#endif

                Rcpp::IntegerVector primary_extract(total_len);
                Index_ current = 0;
                for (const auto& p : to_populate) {
                    Index_ chunk_start = my_chunk_ticks[p.first];
                    Index_ chunk_len = my_chunk_ticks[p.first + 1] - chunk_start;
                    auto start = primary_extract.begin() + current;
                    std::iota(start, start + chunk_len, chunk_start + 1);
                    current += chunk_len;
                }

                my_extract_args[static_cast<int>(!my_row)] = primary_extract;
                auto obj = my_sparse_extractor(my_matrix, my_extract_args);
                parse_sparse_matrix(obj, my_row, my_chunk_value_ptrs, my_chunk_index_ptrs, my_chunk_numbers.data());

                current = 0;
                for (const auto& p : to_populate) {
                    Index_ chunk_len = my_chunk_ticks[p.first + 1] - my_chunk_ticks[p.first];
                    std::copy_n(my_chunk_numbers.begin() + current, chunk_len, p.second->number);
                    current += chunk_len;
                }

#ifdef TATAMI_R_PARALLELIZE_UNKNOWN 
                });
#endif
            }
        );
    }
};

template<bool solo_, bool oracle_, typename Index_, typename CachedValue_, typename CachedIndex_>
using SparseCore = typename std::conditional<solo_,
    SoloSparseCore<oracle_, Index_, CachedValue_, CachedIndex_>,
    typename std::conditional<oracle_,
        OracularSparseCore<Index_, CachedValue_, CachedIndex_>,
        MyopicSparseCore<Index_, CachedValue_, CachedIndex_>
    >::type
>::type;

/******************************
 *** Pure sparse extractors ***
 ******************************/

template<bool solo_, bool oracle_, typename Value_, typename Index_, typename CachedValue_, typename CachedIndex_>
class SparseFull : public tatami::SparseExtractor<oracle_, Value_, Index_> {
public:
    SparseFull(
        const Rcpp::RObject& matrix, 
        const Rcpp::Function& sparse_extractor,
        bool row,
        tatami::MaybeOracle<oracle_, Index_> oracle,
        Index_ non_target_dim,
        Index_ max_target_chunk_length, 
        const std::vector<Index_>& ticks,
        const std::vector<Index_>& map,
        const tatami_chunked::SlabCacheStats& stats,
        bool needs_value,
        bool needs_index) : 
        my_core(
            matrix,
            sparse_extractor,
            row,
            std::move(oracle),
            [&]() {
                Rcpp::IntegerVector output(non_target_dim);
                std::iota(output.begin(), output.end(), 1);
                return output;
            }(),
            max_target_chunk_length,
            ticks,
            map,
            stats,
            needs_value,
            needs_index
        ),
        my_non_target_dim(non_target_dim),
        my_needs_value(needs_value),
        my_needs_index(needs_index)
    {}

private:
    SparseCore<solo_, oracle_, Index_, CachedValue_, CachedIndex_> my_core;
    Index_ my_non_target_dim;
    bool my_needs_value, my_needs_index;

public:
    tatami::SparseRange<Value_, Index_> fetch(Index_ i, Value_* value_buffer, Index_* index_buffer) {
        auto res = my_core.fetch_raw(i);
        const auto& slab = *(res.first);
        Index_ offset = res.second;

        tatami::SparseRange<Value_, Index_> output(slab.number[offset]);
        if (my_needs_value) {
            std::copy_n(slab.values[offset], my_non_target_dim, value_buffer);
            output.value = value_buffer;
        }

        if (my_needs_index) {
            std::copy_n(slab.indices[offset], my_non_target_dim, index_buffer);
            output.index = index_buffer;
        }

        return output;
    }
};

template<bool solo_, bool oracle_, typename Value_, typename Index_, typename CachedValue_, typename CachedIndex_>
class SparseBlock : public tatami::SparseExtractor<oracle_, Value_, Index_> {
public:
    SparseBlock(
        const Rcpp::RObject& matrix, 
        const Rcpp::Function& sparse_extractor,
        bool row,
        tatami::MaybeOracle<oracle_, Index_> oracle,
        Index_ block_start,
        Index_ block_length,
        Index_ max_target_chunk_length, 
        const std::vector<Index_>& ticks,
        const std::vector<Index_>& map,
        const tatami_chunked::SlabCacheStats& stats,
        bool needs_value,
        bool needs_index) : 
        my_core(
            matrix,
            sparse_extractor,
            row,
            std::move(oracle),
            [&]() {
                Rcpp::IntegerVector output(block_length);
                std::iota(output.begin(), output.end(), block_start + 1);
                return output;
            }(),
            max_target_chunk_length,
            ticks,
            map,
            stats,
            needs_value,
            needs_index
        ),
        my_block_start(block_start),
        my_needs_value(needs_value),
        my_needs_index(needs_index)
    {}

private:
    SparseCore<solo_, oracle_, Index_, CachedValue_, CachedIndex_> my_core;
    Index_ my_block_start; 
    bool my_needs_value, my_needs_index;

public:
    tatami::SparseRange<Value_, Index_> fetch(Index_ i, Value_* value_buffer, Index_* index_buffer) {
        auto res = my_core.fetch_raw(i);
        const auto& slab = *(res.first);
        Index_ offset = res.second;

        tatami::SparseRange<Value_, Index_> output(slab.number[offset]);
        if (my_needs_value) {
            std::copy_n(slab.values[offset], output.number, value_buffer); 
            output.value = value_buffer;
        }

        if (my_needs_index) {
            auto iptr = slab.indices[offset];
            for (Index_ i = 0; i < output.number; ++i) {
                index_buffer[i] = static_cast<Index_>(iptr[i]) + my_block_start;
            }
            output.index = index_buffer;
        }

        return output;
    }
};

template<bool solo_, bool oracle_, typename Value_, typename Index_, typename CachedValue_, typename CachedIndex_>
class SparseIndexed : public tatami::SparseExtractor<oracle_, Value_, Index_> {
public:
    SparseIndexed(
        const Rcpp::RObject& matrix, 
        const Rcpp::Function& sparse_extractor,
        bool row,
        tatami::MaybeOracle<oracle_, Index_> oracle,
        tatami::VectorPtr<Index_> indices_ptr,
        Index_ max_target_chunk_length, 
        const std::vector<Index_>& ticks,
        const std::vector<Index_>& map,
        const tatami_chunked::SlabCacheStats& stats,
        bool needs_value,
        bool needs_index) : 
        my_core(
            matrix,
            sparse_extractor,
            row,
            std::move(oracle),
            [&]() {
                Rcpp::IntegerVector output(indices_ptr->begin(), indices_ptr->end());
                for (auto& i : output) {
                    ++i;
                }
                return output;
            }(),
            max_target_chunk_length,
            ticks,
            map,
            stats,
            needs_value,
            needs_index
        ),
        my_indices_ptr(std::move(indices_ptr)),
        my_needs_value(needs_value),
        my_needs_index(needs_index)
    {}

private:
    SparseCore<solo_, oracle_, Index_, CachedValue_, CachedIndex_> my_core;
    tatami::VectorPtr<Index_> my_indices_ptr;
    bool my_needs_value, my_needs_index;

public:
    tatami::SparseRange<Value_, Index_> fetch(Index_ i, Value_* value_buffer, Index_* index_buffer) {
        auto res = my_core.fetch_raw(i);
        const auto& slab = *(res.first);
        Index_ offset = res.second;

        tatami::SparseRange<Value_, Index_> output(slab.number[offset]);
        if (my_needs_value) {
            std::copy_n(slab.values[offset], output.number, value_buffer); 
            output.value = value_buffer;
        }

        if (my_needs_index) {
            auto iptr = slab.indices[offset];
            const auto& indices = *my_indices_ptr;
            for (Index_ i = 0; i < output.number; ++i) {
                index_buffer[i] = indices[iptr[i]];
            }
            output.index = index_buffer;
        }

        return output;
    }
};

/***********************************
 *** Densified sparse extractors ***
 ***********************************/

template<typename Slab_, typename Value_, typename Index_>
const Value_* densify(const Slab_& slab, Index_ offset, size_t non_target_length, Value_* buffer) {
    auto vptr = slab.values[offset];
    auto iptr = slab.indices[offset];
    std::fill_n(buffer, non_target_length, 0);
    for (Index_ i = 0, end = slab.number[offset]; i < end; ++i, ++vptr, ++iptr) {
        buffer[*iptr] = *vptr;
    }
    return buffer;
}

template<bool solo_, bool oracle_, typename Value_, typename Index_, typename CachedValue_, typename CachedIndex_>
class DensifiedSparseFull : public tatami::DenseExtractor<oracle_, Value_, Index_> {
public:
    DensifiedSparseFull(
        const Rcpp::RObject& matrix, 
        const Rcpp::Function& sparse_extractor,
        bool row,
        tatami::MaybeOracle<oracle_, Index_> oracle,
        Index_ non_target_dim,
        Index_ max_target_chunk_length, 
        const std::vector<Index_>& ticks,
        const std::vector<Index_>& map,
        const tatami_chunked::SlabCacheStats& stats) :
        my_core(
            matrix,
            sparse_extractor,
            row,
            std::move(oracle),
            [&]() {
                Rcpp::IntegerVector output(non_target_dim);
                std::iota(output.begin(), output.end(), 1);
                return output;
            }(),
            max_target_chunk_length,
            ticks,
            map,
            stats,
            true,
            true
        ),
        my_non_target_dim(non_target_dim)
    {}

private:
    SparseCore<solo_, oracle_, Index_, CachedValue_, CachedIndex_> my_core;
    size_t my_non_target_dim;

public:
    const Value_* fetch(Index_ i, Value_* buffer) {
        auto res = my_core.fetch_raw(i);
        return densify(*(res.first), res.second, my_non_target_dim, buffer);
    }
};

template<bool solo_, bool oracle_, typename Value_, typename Index_, typename CachedValue_, typename CachedIndex_>
class DensifiedSparseBlock : public tatami::DenseExtractor<oracle_, Value_, Index_> {
public:
    DensifiedSparseBlock(
        const Rcpp::RObject& matrix, 
        const Rcpp::Function& sparse_extractor,
        bool row,
        tatami::MaybeOracle<oracle_, Index_> oracle,
        Index_ block_start,
        Index_ block_length,
        Index_ max_target_chunk_length, 
        const std::vector<Index_>& ticks,
        const std::vector<Index_>& map,
        const tatami_chunked::SlabCacheStats& stats) :
        my_core(
            matrix,
            sparse_extractor,
            row,
            std::move(oracle),
            [&]() {
                Rcpp::IntegerVector output(block_length);
                std::iota(output.begin(), output.end(), block_start + 1);
                return output;
            }(),
            max_target_chunk_length,
            ticks,
            map,
            stats,
            true,
            true
        ),
        my_block_length(block_length)
    {}

private:
    SparseCore<solo_, oracle_, Index_, CachedValue_, CachedIndex_> my_core;
    size_t my_block_length;

public:
    const Value_* fetch(Index_ i, Value_* buffer) {
        auto res = my_core.fetch_raw(i);
        return densify(*(res.first), res.second, my_block_length, buffer);
    }
};

template<bool solo_, bool oracle_, typename Value_, typename Index_, typename CachedValue_, typename CachedIndex_>
class DensifiedSparseIndexed : public tatami::DenseExtractor<oracle_, Value_, Index_> {
public:
    DensifiedSparseIndexed(
        const Rcpp::RObject& matrix, 
        const Rcpp::Function& sparse_extractor,
        bool row,
        tatami::MaybeOracle<oracle_, Index_> oracle,
        tatami::VectorPtr<Index_> idx_ptr,
        Index_ max_target_chunk_length, 
        const std::vector<Index_>& ticks,
        const std::vector<Index_>& map,
        const tatami_chunked::SlabCacheStats& stats) :
        my_core( 
            matrix,
            sparse_extractor,
            row,
            std::move(oracle),
            [&]() {
                Rcpp::IntegerVector output(idx_ptr->begin(), idx_ptr->end());
                for (auto& i : output) {
                    ++i;
                }
                return output;
            }(),
            max_target_chunk_length,
            ticks,
            map,
            stats,
            true,
            true
        ),
        my_num_indices(idx_ptr->size())
    {}

private:
    SparseCore<solo_, oracle_, Index_, CachedValue_, CachedIndex_> my_core;
    size_t my_num_indices;

public:
    const Value_* fetch(Index_ i, Value_* buffer) {
        auto res = my_core.fetch_raw(i);
        return densify(*(res.first), res.second, my_num_indices, buffer);
    }
};

}

}

#endif
