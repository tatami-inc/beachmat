#ifndef TATAMI_CHUNKED_CUSTOM_CHUNK_COORDINATOR_HPP
#define TATAMI_CHUNKED_CUSTOM_CHUNK_COORDINATOR_HPP

#include "tatami/tatami.hpp"
#include "DenseSlabFactory.hpp"
#include "SparseSlabFactory.hpp"
#include "OracularSlabCache.hpp"
#include "OracularSubsettedSlabCache.hpp"
#include "ChunkDimensionStats.hpp"

#include <vector>

namespace tatami_chunked {

namespace CustomChunkedMatrix_internal {

/******************
 *** Workspaces ***
 ******************/

template<typename CachedValue_>
using DenseSingleWorkspace = std::vector<CachedValue_>;

template<typename CachedValue_, typename Index_>
class SparseSingleWorkspace {
public:
    SparseSingleWorkspace(size_t primary_chunkdim, size_t secondary_chunkdim, bool needs_value, bool needs_index) : my_number(primary_chunkdim) {
        size_t total_size = primary_chunkdim * secondary_chunkdim;
        if (needs_value) {
            my_value_pool.resize(total_size);
            my_values.reserve(primary_chunkdim);
            auto vptr = my_value_pool.data();
            for (size_t p = 0; p < primary_chunkdim; ++p, vptr += secondary_chunkdim) {
                my_values.push_back(vptr);
            }
        }
        if (needs_index) {
            my_index_pool.resize(total_size);
            my_indices.reserve(primary_chunkdim);
            auto iptr = my_index_pool.data();
            for (size_t p = 0; p < primary_chunkdim; ++p, iptr += secondary_chunkdim) {
                my_indices.push_back(iptr);
            }
        }
    }

    // Delete the copy constructors as we're passing out pointers.
    SparseSingleWorkspace(const SparseSingleWorkspace&) = delete;
    SparseSingleWorkspace& operator=(const SparseSingleWorkspace&) = delete;

    // Move constructors are okay though.
    SparseSingleWorkspace(SparseSingleWorkspace&&) = default;
    SparseSingleWorkspace& operator=(SparseSingleWorkspace&&) = default;

private:
    std::vector<CachedValue_> my_value_pool;
    std::vector<Index_> my_index_pool;
    std::vector<CachedValue_*> my_values;
    std::vector<Index_*> my_indices;
    std::vector<Index_> my_number;

public:
    std::vector<CachedValue_*>& get_values() {
        return my_values;
    }

    std::vector<Index_*>& get_indices() {
        return my_indices;
    }

    std::vector<Index_>& get_number() {
        return my_number;
    }
};

/*******************
 *** Coordinator ***
 *******************/

template<typename Index_, bool sparse_, class Chunk_>
class ChunkCoordinator {
public:
    ChunkCoordinator(ChunkDimensionStats<Index_> row_stats, ChunkDimensionStats<Index_> col_stats, std::vector<Chunk_> chunk_array, bool row_major) :
        my_row_stats(std::move(row_stats)), my_col_stats(std::move(col_stats)), my_chunk_array(std::move(chunk_array)), my_row_major(row_major)
    {
        if (static_cast<size_t>(my_row_stats.num_chunks) * static_cast<size_t>(my_col_stats.num_chunks) != my_chunk_array.size()) {
            throw std::runtime_error("length of 'chunks' should be equal to the product of the number of chunks along each row and column");
        }
    }

private:
    ChunkDimensionStats<Index_> my_row_stats;
    ChunkDimensionStats<Index_> my_col_stats;
    std::vector<Chunk_> my_chunk_array;
    bool my_row_major;

public:
    // Number of chunks along the rows is equal to the number of chunks for
    // each column, and vice versa; hence the flipped definitions.
    Index_ get_num_chunks_per_row() const {
        return my_col_stats.num_chunks;
    }

    Index_ get_num_chunks_per_column() const {
        return my_row_stats.num_chunks;
    }

    Index_ get_nrow() const {
        return my_row_stats.dimension_extent;
    }

    Index_ get_ncol() const {
        return my_col_stats.dimension_extent;
    }

    bool prefer_rows_internal() const {
        // Prefer rows if we have to extract fewer chunks per row.
        return get_num_chunks_per_column() > get_num_chunks_per_row(); 
    }

    Index_ get_chunk_nrow() const {
        return my_row_stats.chunk_length;
    }

    Index_ get_chunk_ncol() const {
        return my_col_stats.chunk_length;
    }

public:
    Index_ get_secondary_dim(bool row) const {
        if (row) {
            return my_col_stats.dimension_extent;
        } else {
            return my_row_stats.dimension_extent;
        }
    }

    Index_ get_primary_chunkdim(bool row) const {
        if (row) {
            return my_row_stats.chunk_length;
        } else {
            return my_col_stats.chunk_length;
        }
    }

    Index_ get_secondary_chunkdim(bool row) const {
        if (row) {
            return my_col_stats.chunk_length;
        } else {
            return my_row_stats.chunk_length;
        }
    }

    // Overload that handles the truncated chunk at the bottom/right edges of each matrix.
    Index_ get_primary_chunkdim(bool row, Index_ chunk_id) const {
        return get_chunk_length(row ? my_row_stats : my_col_stats, chunk_id);
    }

private:
    std::pair<size_t, size_t> offset_and_increment(bool row, Index_ chunk_id) const {
        size_t num_chunks = (my_row_major ? get_num_chunks_per_row() : get_num_chunks_per_column()); // use size_t to avoid overflow.
        if (row == my_row_major) {
            return std::pair<size_t, size_t>(static_cast<size_t>(chunk_id) * num_chunks, 1);
        } else {
            return std::pair<size_t, size_t>(chunk_id, num_chunks);
        }
    }

    template<class ExtractFunction_>
    void extract_secondary_block(
        bool row,
        Index_ chunk_id, 
        Index_ secondary_block_start, 
        Index_ secondary_block_length, 
        ExtractFunction_ extract)
    const {
        auto secondary_chunkdim = get_secondary_chunkdim(row);
        Index_ start_chunk_index = secondary_block_start / secondary_chunkdim;
        Index_ secondary_start_pos = start_chunk_index * secondary_chunkdim;
        Index_ secondary_block_end = secondary_block_start + secondary_block_length;
        Index_ end_chunk_index = secondary_block_end / secondary_chunkdim + (secondary_block_end % secondary_chunkdim > 0); // i.e., integer ceiling.

        auto oi = offset_and_increment(row, chunk_id);
        auto offset = std::get<0>(oi);
        auto increment = std::get<1>(oi);
        offset += increment * static_cast<size_t>(start_chunk_index); // size_t to avoid integer overflow.

        for (Index_ c = start_chunk_index; c < end_chunk_index; ++c) {
            const auto& chunk = my_chunk_array[offset];
            Index_ from = (c == start_chunk_index ? secondary_block_start - secondary_start_pos : 0);
            Index_ to = (c + 1 == end_chunk_index ? secondary_block_end - secondary_start_pos : secondary_chunkdim);
            Index_ len = to - from;

            // No need to protect against a zero length, as it should be impossible
            // here (otherwise, start_chunk_index == end_chunk_index and we'd never iterate).
            if constexpr(sparse_) {
                extract(chunk, from, len, secondary_start_pos);
            } else {
                extract(chunk, from, len);
            }

            secondary_start_pos += to; // yes, this is deliberate; '+ to' means that either we add 'secondary_chunkdim' or set it to 'secondary_block_end', the latter of which avoids overflow.
            offset += increment;
        }
    }

    template<class ExtractFunction_>
    void extract_secondary_index(
        bool row,
        Index_ chunk_id, 
        const std::vector<Index_>& secondary_indices, 
        std::vector<Index_>& chunk_indices_buffer, 
        ExtractFunction_ extract)
    const {
        if (secondary_indices.empty()) {
            return;
        }

        auto secondary_chunkdim = get_secondary_chunkdim(row);
        Index_ start_chunk_index = secondary_indices.front() / secondary_chunkdim; // 'secondary_indices' is guaranteed to be non-empty at this point.
        Index_ secondary_start_pos = start_chunk_index * secondary_chunkdim;

        auto oi = offset_and_increment(row, chunk_id);
        auto offset = std::get<0>(oi);
        auto increment = std::get<1>(oi);
        offset += increment * static_cast<size_t>(start_chunk_index); // size_t to avoid integer overflow.

        auto secondary_dim = get_secondary_dim(row);
        auto iIt = secondary_indices.begin();
        auto iEnd = secondary_indices.end();
        while (iIt != iEnd) {
            const auto& chunk = my_chunk_array[offset];

            Index_ secondary_end_pos = std::min(secondary_dim - secondary_start_pos, secondary_chunkdim) + secondary_start_pos; // this convoluted method avoids overflow.
            chunk_indices_buffer.clear();
            while (iIt != iEnd && *iIt < secondary_end_pos) {
                chunk_indices_buffer.push_back(*iIt - secondary_start_pos);
                ++iIt;
            }

            if (!chunk_indices_buffer.empty()) {
                if constexpr(sparse_) {
                    extract(chunk, chunk_indices_buffer, secondary_start_pos);
                } else {
                    extract(chunk, chunk_indices_buffer);
                }
            }

            secondary_start_pos = secondary_end_pos;
            offset += increment;
        }
    }

    typedef typename Chunk_::Workspace ChunkWork;
    typedef typename Chunk_::value_type ChunkValue;
    typedef typename std::conditional<sparse_, typename SparseSlabFactory<ChunkValue, Index_>::Slab, typename DenseSlabFactory<ChunkValue>::Slab>::type Slab;
    typedef typename std::conditional<sparse_, SparseSingleWorkspace<ChunkValue, Index_>, DenseSingleWorkspace<ChunkValue> >::type SingleWorkspace;

public:
    // Extract a single element of the primary dimension, using a contiguous
    // block on the secondary dimension. 
    //
    // Unfortunately, we can't just re-use the fetch_block() functions with a
    // length of 1, because some chunks do not support partial extraction; this
    // requires special handling to extract the full chunk into 'tmp_work', and
    // then pull out what we need into 'final_slab'.
    //
    // Even the use_subset=true chunks that do support partial extraction
    // require a workspace involving the full chunk size, so we end up needing
    // the 'tmp_work' workspace anyway.
    std::pair<const Slab*, Index_> fetch_single(
        bool row,
        Index_ i,
        Index_ secondary_block_start, 
        Index_ secondary_block_length, 
        ChunkWork& chunk_workspace,
        SingleWorkspace& tmp_work,
        Slab& final_slab)
    const {
        Index_ primary_chunkdim = get_primary_chunkdim(row);
        Index_ chunk_id = i / primary_chunkdim;
        Index_ chunk_offset = i % primary_chunkdim;

        if constexpr(sparse_) {
            auto& final_num = *final_slab.number;
            final_num = 0;
            bool needs_value = !final_slab.values.empty();
            bool needs_index = !final_slab.indices.empty();

            extract_secondary_block(
                row, chunk_id, secondary_block_start, secondary_block_length, 
                [&](const Chunk_& chunk, Index_ from, Index_ len, Index_ secondary_start_pos) {
                    auto& tmp_values = tmp_work.get_values();
                    auto& tmp_indices = tmp_work.get_indices();
                    auto& tmp_number = tmp_work.get_number();
                    std::fill_n(tmp_number.begin(), primary_chunkdim, 0);

                    if constexpr(Chunk_::use_subset) {
                        chunk.extract(row, chunk_offset, 1, from, len, chunk_workspace, tmp_values, tmp_indices, tmp_number.data(), secondary_start_pos);
                    } else {
                        chunk.extract(row, from, len, chunk_workspace, tmp_values, tmp_indices, tmp_number.data(), secondary_start_pos);
                    }

                    auto count = tmp_number[chunk_offset];
                    if (needs_value) {
                        std::copy_n(tmp_values[chunk_offset], count, final_slab.values[0] + final_num);
                    }
                    if (needs_index) {
                        std::copy_n(tmp_indices[chunk_offset], count, final_slab.indices[0] + final_num);
                    }
                    final_num += count;
                }
            );

        } else {
            auto final_slab_ptr = final_slab.data;
            auto tmp_buffer_ptr = tmp_work.data();

            extract_secondary_block(
                row, chunk_id, secondary_block_start, secondary_block_length, 
                [&](const Chunk_& chunk, Index_ from, Index_ len) {
                    if constexpr(Chunk_::use_subset) {
                        chunk.extract(row, chunk_offset, 1, from, len, chunk_workspace, tmp_buffer_ptr, len);
                    } else {
                        chunk.extract(row, from, len, chunk_workspace, tmp_buffer_ptr, len);
                    }

                    size_t tmp_offset = static_cast<size_t>(len) * static_cast<size_t>(chunk_offset);
                    std::copy_n(tmp_buffer_ptr + tmp_offset, len, final_slab_ptr);
                    final_slab_ptr += len;
                }
            );
        }

        return std::make_pair(&final_slab, static_cast<Index_>(0));
    }

    // Extract a single element of the primary dimension, using an indexed
    // subset on the secondary dimension.
    std::pair<const Slab*, Index_> fetch_single(
        bool row,
        Index_ i,
        const std::vector<Index_>& secondary_indices, 
        std::vector<Index_>& chunk_indices_buffer,
        ChunkWork& chunk_workspace,
        SingleWorkspace& tmp_work,
        Slab& final_slab)
    const {
        Index_ primary_chunkdim = get_primary_chunkdim(row);
        Index_ chunk_id = i / primary_chunkdim;
        Index_ chunk_offset = i % primary_chunkdim;

        if constexpr(sparse_) {
            auto& final_num = *final_slab.number;
            final_num = 0;
            bool needs_value = !final_slab.values.empty();
            bool needs_index = !final_slab.indices.empty();

            extract_secondary_index(
                row, chunk_id, secondary_indices, chunk_indices_buffer,
                [&](const Chunk_& chunk, const std::vector<Index_>& chunk_indices, Index_ secondary_start_pos) {
                    auto& tmp_values = tmp_work.get_values();
                    auto& tmp_indices = tmp_work.get_indices();
                    auto& tmp_number = tmp_work.get_number();
                    std::fill_n(tmp_number.begin(), primary_chunkdim, 0);

                    if constexpr(Chunk_::use_subset) {
                        chunk.extract(row, chunk_offset, 1, chunk_indices, chunk_workspace, tmp_values, tmp_indices, tmp_number.data(), secondary_start_pos);
                    } else {
                        chunk.extract(row, chunk_indices, chunk_workspace, tmp_values, tmp_indices, tmp_number.data(), secondary_start_pos);
                    }

                    auto count = tmp_number[chunk_offset];
                    if (needs_value) {
                        std::copy_n(tmp_values[chunk_offset], count, final_slab.values[0] + final_num);
                    }
                    if (needs_index) {
                        std::copy_n(tmp_indices[chunk_offset], count, final_slab.indices[0] + final_num);
                    }
                    final_num += count;
                }
            );

        } else {
            auto final_slab_ptr = final_slab.data;
            auto tmp_buffer_ptr = tmp_work.data();

            extract_secondary_index(
                row, chunk_id, secondary_indices, chunk_indices_buffer,
                [&](const Chunk_& chunk, const std::vector<Index_>& chunk_indices) {
                    size_t nidx = chunk_indices.size();
                    if constexpr(Chunk_::use_subset) {
                        chunk.extract(row, chunk_offset, 1, chunk_indices, chunk_workspace, tmp_buffer_ptr, nidx);
                    } else {
                        chunk.extract(row, chunk_indices, chunk_workspace, tmp_buffer_ptr, nidx);
                    }

                    size_t tmp_offset = nidx * static_cast<size_t>(chunk_offset);
                    std::copy_n(tmp_buffer_ptr + tmp_offset, nidx, final_slab_ptr);
                    final_slab_ptr += nidx;
                }
            );
        }

        return std::make_pair(&final_slab, static_cast<Index_>(0));
    }

private:
    // Extract a contiguous block of the primary dimension, using a contiguous block on the secondary dimension.
    void fetch_block(
        bool row,
        Index_ chunk_id, 
        Index_ chunk_offset, 
        Index_ chunk_length, 
        Index_ secondary_block_start, 
        Index_ secondary_block_length, 
        Slab& slab, 
        ChunkWork& chunk_workspace) 
    const {
        if constexpr(sparse_) {
            std::fill_n(slab.number, get_primary_chunkdim(row), 0);

            extract_secondary_block(
                row, chunk_id, secondary_block_start, secondary_block_length, 
                [&](const Chunk_& chunk, Index_ from, Index_ len, Index_ secondary_start_pos) {
                    if constexpr(Chunk_::use_subset) {
                        chunk.extract(row, chunk_offset, chunk_length, from, len, chunk_workspace, slab.values, slab.indices, slab.number, secondary_start_pos);
                    } else {
                        chunk.extract(row, from, len, chunk_workspace, slab.values, slab.indices, slab.number, secondary_start_pos);
                    }
                }
            );

        } else {
            auto slab_ptr = slab.data;
            size_t stride = secondary_block_length;

            extract_secondary_block(
                row, chunk_id, secondary_block_start, secondary_block_length, 
                [&](const Chunk_& chunk, Index_ from, Index_ len) {
                    if constexpr(Chunk_::use_subset) {
                        chunk.extract(row, chunk_offset, chunk_length, from, len, chunk_workspace, slab_ptr, stride);
                    } else {
                        chunk.extract(row, from, len, chunk_workspace, slab_ptr, stride);
                    }
                    slab_ptr += len;
                }
            );
        }
    }

    // Extract a contiguous block of the primary dimension, using an indexed subset on the secondary dimension.
    void fetch_block(
        bool row,
        Index_ chunk_id, 
        Index_ chunk_offset, 
        Index_ chunk_length, 
        const std::vector<Index_>& secondary_indices, 
        std::vector<Index_>& chunk_indices_buffer,
        Slab& slab, 
        ChunkWork& chunk_workspace) 
    const {
        if constexpr(sparse_) {
            std::fill_n(slab.number, get_primary_chunkdim(row), 0);

            extract_secondary_index(
                row, chunk_id, secondary_indices, chunk_indices_buffer,
                [&](const Chunk_& chunk, const std::vector<Index_>& chunk_indices, Index_ secondary_start_pos) {
                    if constexpr(Chunk_::use_subset) {
                        chunk.extract(row, chunk_offset, chunk_length, chunk_indices, chunk_workspace, slab.values, slab.indices, slab.number, secondary_start_pos);
                    } else {
                        chunk.extract(row, chunk_indices, chunk_workspace, slab.values, slab.indices, slab.number, secondary_start_pos);
                    }
                }
            );

        } else {
            auto slab_ptr = slab.data;
            size_t stride = secondary_indices.size();

            extract_secondary_index(
                row, chunk_id, secondary_indices, chunk_indices_buffer,
                [&](const Chunk_& chunk, const std::vector<Index_>& chunk_indices) {
                    if constexpr(Chunk_::use_subset) {
                        chunk.extract(row, chunk_offset, chunk_length, chunk_indices, chunk_workspace, slab_ptr, stride);
                    } else {
                        chunk.extract(row, chunk_indices, chunk_workspace, slab_ptr, stride);
                    }
                    slab_ptr += chunk_indices.size();
                }
            );
        }
    }

private:
    // Extract an indexed subset of the primary dimension, using a contiguous block on the secondary dimension.
    void fetch_index(
        bool row,
        Index_ chunk_id,
        const std::vector<Index_>& primary_indices, 
        Index_ secondary_block_start, 
        Index_ secondary_block_length, 
        Slab& slab, 
        ChunkWork& chunk_workspace) 
    const {
        if constexpr(sparse_) {
            std::fill_n(slab.number, get_primary_chunkdim(row), 0);

            extract_secondary_block(
                row, chunk_id, secondary_block_start, secondary_block_length,
                [&](const Chunk_& chunk, Index_ from, Index_ len, Index_ secondary_start_pos) {
                    if constexpr(Chunk_::use_subset) {
                        chunk.extract(row, primary_indices, from, len, chunk_workspace, slab.values, slab.indices, slab.number, secondary_start_pos);
                    } else {
                        chunk.extract(row, from, len, chunk_workspace, slab.values, slab.indices, slab.number, secondary_start_pos);
                    }
                }
            );

        } else {
            auto slab_ptr = slab.data;
            size_t stride = secondary_block_length;

            extract_secondary_block(
                row, chunk_id, secondary_block_start, secondary_block_length,
                [&](const Chunk_& chunk, Index_ from, Index_ len) {
                    if constexpr(Chunk_::use_subset) {
                        chunk.extract(row, primary_indices, from, len, chunk_workspace, slab_ptr, stride);
                    } else {
                        chunk.extract(row, from, len, chunk_workspace, slab_ptr, stride);
                    }
                    slab_ptr += len;
                }
            );
        }
    }

    // Extract an indexed subset of the primary dimension, using an indexed subset on the secondary dimension.
    void fetch_index(
        bool row,
        Index_ chunk_id,
        const std::vector<Index_>& primary_indices,
        const std::vector<Index_>& secondary_indices,
        std::vector<Index_>& chunk_indices_buffer,
        Slab& slab, 
        ChunkWork& chunk_workspace) 
    const {
        if constexpr(sparse_) {
            std::fill_n(slab.number, get_primary_chunkdim(row), 0);

            extract_secondary_index(
                row, chunk_id, secondary_indices, chunk_indices_buffer,
                [&](const Chunk_& chunk, const std::vector<Index_>& chunk_indices, Index_ secondary_start_pos) {
                    if constexpr(Chunk_::use_subset) {
                        chunk.extract(row, primary_indices, chunk_indices, chunk_workspace, slab.values, slab.indices, slab.number, secondary_start_pos);
                    } else {
                        chunk.extract(row, chunk_indices, chunk_workspace, slab.values, slab.indices, slab.number, secondary_start_pos);
                    }
                }
            );

        } else {
            auto slab_ptr = slab.data;
            size_t stride = secondary_indices.size();

            extract_secondary_index(
                row, chunk_id, secondary_indices, chunk_indices_buffer, 
                [&](const Chunk_& chunk, const std::vector<Index_>& chunk_indices) {
                    if constexpr(Chunk_::use_subset) {
                        chunk.extract(row, primary_indices, chunk_indices, chunk_workspace, slab_ptr, stride);
                    } else {
                        chunk.extract(row, chunk_indices, chunk_workspace, slab_ptr, stride);
                    }
                    slab_ptr += chunk_indices.size();
                }
            );
        }
    }

public:
    // Obtain the slab containing the 'i'-th element of the primary dimension.
    template<class Cache_, class Factory_>
    std::pair<const Slab*, Index_> fetch_myopic(
        bool row,
        Index_ i, 
        Index_ block_start,
        Index_ block_length,
        ChunkWork& chunk_workspace,
        Cache_& cache,
        Factory_& factory)
    const {
        Index_ primary_chunkdim = get_primary_chunkdim(row);
        Index_ chunk_id = i / primary_chunkdim;
        Index_ chunk_offset = i % primary_chunkdim;
        auto& out = cache.find(
            chunk_id,
            /* create = */ [&]() -> Slab {
                return factory.create();
            },
            /* populate = */ [&](Index_ id, Slab& slab) -> void {
                fetch_block(row, id, 0, get_primary_chunkdim(row, id), block_start, block_length, slab, chunk_workspace);
            }
        );
        return std::make_pair(&out, chunk_offset);
    }

    template<class Cache_, class Factory_>
    std::pair<const Slab*, Index_> fetch_myopic(
        bool row,
        Index_ i, 
        const std::vector<Index_>& indices,
        std::vector<Index_>& tmp_indices,
        ChunkWork& chunk_workspace,
        Cache_& cache,
        Factory_& factory)
    const {
        Index_ primary_chunkdim = get_primary_chunkdim(row);
        Index_ chunk_id = i / primary_chunkdim;
        Index_ chunk_offset = i % primary_chunkdim;
        auto& out = cache.find(
            chunk_id,
            /* create = */ [&]() -> Slab {
                return factory.create();
            },
            /* populate = */ [&](Index_ id, Slab& slab) -> void {
                fetch_block(row, id, 0, get_primary_chunkdim(row, id), indices, tmp_indices, slab, chunk_workspace);
            }
        );
        return std::make_pair(&out, chunk_offset);
    }

private:
    template<class Cache_, class Factory_, typename PopulateBlock_, typename PopulateIndex_>
    std::pair<const Slab*, Index_> fetch_oracular_core(
        bool row,
        Cache_& cache,
        Factory_& factory,
        PopulateBlock_ populate_block,
        PopulateIndex_ populate_index)
    const {
        Index_ primary_chunkdim = get_primary_chunkdim(row);
        if constexpr(Chunk_::use_subset) {
            return cache.next(
                /* identify = */ [&](Index_ i) -> std::pair<Index_, Index_> {
                    return std::pair<Index_, Index_>(i / primary_chunkdim, i % primary_chunkdim);
                },
                /* create = */ [&]() -> Slab {
                    return factory.create();
                },
                /* populate =*/ [&](std::vector<std::tuple<Index_, Slab*, const OracularSubsettedSlabCacheSelectionDetails<Index_>*> >& in_need) -> void {
                    for (const auto& p : in_need) {
                        auto id = std::get<0>(p);
                        auto ptr = std::get<1>(p);
                        auto sub = std::get<2>(p);
                        switch (sub->selection) {
                            case OracularSubsettedSlabCacheSelectionType::FULL:
                                populate_block(id, 0, get_primary_chunkdim(row, id), *ptr);
                                break;
                            case OracularSubsettedSlabCacheSelectionType::BLOCK:
                                populate_block(id, sub->block_start, sub->block_length, *ptr);
                                break;
                            case OracularSubsettedSlabCacheSelectionType::INDEX:
                                populate_index(id, sub->indices, *ptr);
                                break;
                        }
                    }
                }
            );

        } else {
            return cache.next(
                /* identify = */ [&](Index_ i) -> std::pair<Index_, Index_> {
                    return std::pair<Index_, Index_>(i / primary_chunkdim, i % primary_chunkdim);
                },
                /* create = */ [&]() -> Slab {
                    return factory.create();
                },
                /* populate =*/ [&](std::vector<std::pair<Index_, Slab*> >& to_populate) -> void {
                    for (auto& p : to_populate) {
                        populate_block(p.first, 0, get_primary_chunkdim(row, p.first), *(p.second));
                    }
                }
            );
        }
    }

public:
    template<class Cache_, class Factory_>
    std::pair<const Slab*, Index_> fetch_oracular(
        bool row,
        Index_ block_start,
        Index_ block_length,
        ChunkWork& chunk_workspace,
        Cache_& cache,
        Factory_& factory)
    const {
        return fetch_oracular_core(
            row,
            cache,
            factory,
            [&](Index_ pid, Index_ pstart, Index_ plen, Slab& slab) {
                fetch_block(row, pid, pstart, plen, block_start, block_length, slab, chunk_workspace);
            },
            [&](Index_ pid, const std::vector<Index_>& pindices, Slab& slab) {
                fetch_index(row, pid, pindices, block_start, block_length, slab, chunk_workspace);
            }
        );
    }

    template<class Cache_, class Factory_>
    std::pair<const Slab*, Index_> fetch_oracular(
        bool row,
        const std::vector<Index_>& indices,
        std::vector<Index_>& chunk_indices_buffer,
        ChunkWork& chunk_workspace,
        Cache_& cache,
        Factory_& factory)
    const {
        return fetch_oracular_core(
            row,
            cache,
            factory,
            [&](Index_ pid, Index_ pstart, Index_ plen, Slab& slab) {
                fetch_block(row, pid, pstart, plen, indices, chunk_indices_buffer, slab, chunk_workspace);
            },
            [&](Index_ pid, const std::vector<Index_>& pindices, Slab& slab) {
                fetch_index(row, pid, pindices, indices, chunk_indices_buffer, slab, chunk_workspace);
            }
        );
    }
};

}

}

#endif
