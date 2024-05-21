#ifndef TATAMI_CHUNKED_CUSTOM_DENSE_CHUNKED_MATRIX_HPP
#define TATAMI_CHUNKED_CUSTOM_DENSE_CHUNKED_MATRIX_HPP

#include "tatami/tatami.hpp"
#include "custom_internals.hpp"
#include "SlabCacheStats.hpp"
#include "DenseSlabFactory.hpp"
#include "LruSlabCache.hpp"
#include "OracularSlabCache.hpp"
#include "OracularSubsettedSlabCache.hpp"

#include <type_traits>
#include <vector>

/**
 * @file CustomDenseChunkedMatrix.hpp
 * @brief Custom dense chunked matrix.
 */

namespace tatami_chunked {

/**
 * @brief Options for data extraction from a `CustomDenseChunkedMatrix`.
 */
struct CustomDenseChunkedMatrixOptions {
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
 * @cond
 */
namespace CustomChunkedMatrix_internal {

/*********************
 **** Base classes ***
 *********************/

template<bool oracle_, typename Value_, typename Index_, typename Chunk_>
class SoloDenseCore {
private:
    const ChunkCoordinator<Index_, false, Chunk_>& my_coordinator;
    typename Chunk_::Workspace my_chunk_workspace;

    tatami::MaybeOracle<oracle_, Index_> my_oracle;
    typename std::conditional<oracle_, size_t, bool>::type my_counter = 0;

    DenseSlabFactory<typename Chunk_::value_type> my_factory;
    typedef typename decltype(my_factory)::Slab Slab;

    // These two instances are not fully allocated Slabs; rather, tmp_solo just
    // holds the content for a single chunk, while final_solo holds the content
    // across chunks but only for the requested dimension element. Both cases
    // are likely to be much smaller than a full Slab, so we're already more
    // memory-efficient than 'require_minimum_cache = true`. 
    DenseSingleWorkspace<typename Chunk_::value_type> my_tmp_solo;
    Slab my_final_solo;

public:
    SoloDenseCore(
        const ChunkCoordinator<Index_, false, Chunk_>& coordinator, 
        [[maybe_unused]] const SlabCacheStats& slab_stats, // for consistency with the other base classes.
        tatami::MaybeOracle<oracle_, Index_> oracle,
        Index_ secondary_length) :
        my_coordinator(coordinator),
        my_oracle(std::move(oracle)),
        my_factory(secondary_length, 1),
        my_tmp_solo(static_cast<size_t>(my_coordinator.get_chunk_nrow()) * static_cast<size_t>(my_coordinator.get_chunk_ncol())),
        my_final_solo(my_factory.create())
    {}

    template<typename ... Args_>
    std::pair<const Slab*, Index_> fetch_raw(bool row, Index_ i, Args_&& ... args) {
        if constexpr(oracle_) {
            i = my_oracle->get(my_counter++);
        }
        return my_coordinator.fetch_single(row, i, std::forward<Args_>(args)..., my_chunk_workspace, my_tmp_solo, my_final_solo);
    }
};

template<typename Value_, typename Index_, typename Chunk_>
class MyopicDenseCore {
private:
    const ChunkCoordinator<Index_, false, Chunk_>& my_coordinator;
    typename Chunk_::Workspace my_chunk_workspace;

    DenseSlabFactory<typename Chunk_::value_type> my_factory;
    typedef typename decltype(my_factory)::Slab Slab;

    LruSlabCache<Index_, Slab> my_cache;

public:
    MyopicDenseCore(
        const ChunkCoordinator<Index_, false, Chunk_>& coordinator,
        const SlabCacheStats& slab_stats, 
        [[maybe_unused]] tatami::MaybeOracle<false, Index_> ora, // for consistency with the other base classes
        [[maybe_unused]] Index_ secondary_length) :
        my_coordinator(coordinator),
        my_factory(slab_stats),
        my_cache(slab_stats.max_slabs_in_cache)
    {}

    template<typename ... Args_>
    std::pair<const Slab*, Index_> fetch_raw(bool row, Index_ i, Args_&& ... args) {
        return my_coordinator.fetch_myopic(row, i, std::forward<Args_>(args)..., my_chunk_workspace, my_cache, my_factory);
    }
};

template<typename Value_, typename Index_, typename Chunk_>
class OracularDenseCore {
private:
    const ChunkCoordinator<Index_, false, Chunk_>& my_coordinator;
    typename Chunk_::Workspace my_chunk_workspace;

    DenseSlabFactory<typename Chunk_::value_type> my_factory;
    typedef typename decltype(my_factory)::Slab Slab;

    typename std::conditional<Chunk_::use_subset, OracularSubsettedSlabCache<Index_, Index_, Slab>, OracularSlabCache<Index_, Index_, Slab> >::type my_cache;

public:
    OracularDenseCore(
        const ChunkCoordinator<Index_, false, Chunk_>& coordinator,
        const SlabCacheStats& slab_stats,
        tatami::MaybeOracle<true, Index_> oracle, 
        [[maybe_unused]] Index_ secondary_length) :
        my_coordinator(coordinator), 
        my_factory(slab_stats),
        my_cache(std::move(oracle), slab_stats.max_slabs_in_cache)
    {}

    template<typename ... Args_>
    std::pair<const Slab*, Index_> fetch_raw(bool row, [[maybe_unused]] Index_ i, Args_&& ... args) {
        return my_coordinator.fetch_oracular(row, std::forward<Args_>(args)..., my_chunk_workspace, my_cache, my_factory);
    }
};

template<bool solo_, bool oracle_, typename Value_, typename Index_, typename Chunk_>
using DenseCore = typename std::conditional<solo_, 
      SoloDenseCore<oracle_, Value_, Index_, Chunk_>,
      typename std::conditional<oracle_,
          OracularDenseCore<Value_, Index_, Chunk_>,
          MyopicDenseCore<Value_, Index_, Chunk_>
      >::type
>::type;

/***********************
 **** Actual classes ***
 ***********************/

template<class Slab_, typename Index_, typename Value_>
const Value_* process_dense_slab(const std::pair<const Slab_*, Index_>& fetched, Value_* buffer, size_t secondary_length) {
    auto ptr = fetched.first->data + static_cast<size_t>(fetched.second) * secondary_length; // cast to size_t to avoid overflow.
    std::copy_n(ptr, secondary_length, buffer);
    return buffer;
}

template<bool solo_, bool oracle_, typename Value_, typename Index_, typename Chunk_>
class DenseFull : public tatami::DenseExtractor<oracle_, Value_, Index_> {
public:
    DenseFull(
        const ChunkCoordinator<Index_, false, Chunk_>& coordinator, 
        const SlabCacheStats& slab_stats,
        bool row,
        tatami::MaybeOracle<oracle_, Index_> oracle) :
        my_row(row),
        my_secondary_dim(coordinator.get_secondary_dim(row)),
        my_core(
            coordinator,
            slab_stats,
            std::move(oracle),
            my_secondary_dim
        )
    {}

    const Value_* fetch(Index_ i, Value_* buffer) {
        auto fetched = my_core.fetch_raw(my_row, i, 0, my_secondary_dim);
        return process_dense_slab(fetched, buffer, my_secondary_dim);
    }

private:
    bool my_row;
    Index_ my_secondary_dim;
    DenseCore<solo_, oracle_, Value_, Index_, Chunk_> my_core;
};

template<bool solo_, bool oracle_, typename Value_, typename Index_, typename Chunk_>
class DenseBlock : public tatami::DenseExtractor<oracle_, Value_, Index_> {
public:
    DenseBlock(
        const ChunkCoordinator<Index_, false, Chunk_>& coordinator, 
        const SlabCacheStats& slab_stats,
        bool row,
        tatami::MaybeOracle<oracle_, Index_> ora, 
        Index_ block_start, 
        Index_ block_length) :
        my_row(row),
        my_block_start(block_start),
        my_block_length(block_length),
        my_core(
            coordinator, 
            slab_stats,
            std::move(ora),
            block_length
        )
    {}

    const Value_* fetch(Index_ i, Value_* buffer) {
        auto fetched = my_core.fetch_raw(my_row, i, my_block_start, my_block_length);
        return process_dense_slab(fetched, buffer, my_block_length);
    }

private:
    bool my_row;
    Index_ my_block_start, my_block_length;
    DenseCore<solo_, oracle_, Value_, Index_, Chunk_> my_core;
};

template<bool solo_, bool oracle_, typename Value_, typename Index_, typename Chunk_>
class DenseIndex : public tatami::DenseExtractor<oracle_, Value_, Index_> {
public:
    DenseIndex(
        const ChunkCoordinator<Index_, false, Chunk_>& coordinator, 
        const SlabCacheStats& slab_stats,
        bool row,
        tatami::MaybeOracle<oracle_, Index_> oracle,
        tatami::VectorPtr<Index_> indices_ptr) :
        my_row(row),
        my_indices_ptr(std::move(indices_ptr)),
        my_core(
            coordinator, 
            slab_stats,
            std::move(oracle),
            my_indices_ptr->size()
        )
    {}

    const Value_* fetch(Index_ i, Value_* buffer) {
        auto fetched = my_core.fetch_raw(my_row, i, *my_indices_ptr, my_tmp_indices);
        return process_dense_slab(fetched, buffer, my_indices_ptr->size());
    }

private:
    bool my_row;
    tatami::VectorPtr<Index_> my_indices_ptr;
    std::vector<Index_> my_tmp_indices;
    DenseCore<solo_, oracle_, Value_, Index_, Chunk_> my_core;
};

}
/**
 * @endcond
 */

/**
 * @brief Matrix of custom dense chunks.
 *
 * @tparam Value_ Numeric type for the matrix value.
 * @tparam Index_ Integer type for the row/column indices.
 * @tparam Chunk_ Class of the chunk, implementing either the `MockSimpleDenseChunk` or `MockSubsetDenseChunk` interfaces.
 *
 * Implements a `Matrix` subclass where data is contained in dense rectangular chunks.
 * These chunks are typically compressed in some manner to reduce memory usage compared to, e.g., a `tatami::DenseMatrix`.
 * On access, the relevant chunks are decompressed and the desired values are extracted.
 * Each dimension should be divided into chunk boundaries at regular intervals starting from zero;
 * this partitions the matrix according to a regular grid where each grid entry is a single chunk of the same size.
 * The exception is for chunks at the non-zero boundaries of the matrix dimensions, which may be truncated.
 */
template<typename Value_, typename Index_, typename Chunk_>
class CustomDenseChunkedMatrix : public tatami::Matrix<Value_, Index_> {
public:
    /**
     * @param mat_nrow Number of rows in the matrix.
     * @param mat_ncol Number of columns in the matrix.
     * @param chunk_nrow Number of rows in each chunk.
     * @param chunk_ncol Number of columns in each chunk.
     * @param chunks Vector containing a two-dimensional array of chunks that cover the entire matrix.
     * This should have length equal to the product of the number of chunks along the rows and columns of the matrix, i.e., `ceil(mat_nrow / chunk_nrow) * ceil(mat_ncol / chunk_ncol)`.
     * @param row_major Whether `chunks` is in row-major format.
     * @param opt Further options for chunked extraction.
     */
    CustomDenseChunkedMatrix(Index_ mat_nrow, Index_ mat_ncol, Index_ chunk_nrow, Index_ chunk_ncol, std::vector<Chunk_> chunks, bool row_major, const CustomDenseChunkedMatrixOptions& opt) : 
        my_coordinator(ChunkDimensionStats<Index_>(mat_nrow, chunk_nrow), ChunkDimensionStats<Index_>(mat_ncol, chunk_ncol), std::move(chunks), row_major),
        my_cache_size_in_elements(opt.maximum_cache_size / sizeof(typename Chunk_::value_type)),
        my_require_minimum_cache(opt.require_minimum_cache)
    {}

private:
    CustomChunkedMatrix_internal::ChunkCoordinator<Index_, false, Chunk_> my_coordinator;
    size_t my_cache_size_in_elements;
    bool my_require_minimum_cache;

public:
    Index_ nrow() const { 
        return my_coordinator.get_nrow(); 
    }

    Index_ ncol() const { 
        return my_coordinator.get_ncol(); 
    }

    bool prefer_rows() const { 
        return my_coordinator.prefer_rows_internal();
    }

    bool uses_oracle(bool) const { 
        return true; 
    }

    double prefer_rows_proportion() const { 
        return static_cast<double>(my_coordinator.prefer_rows_internal());
    }

    bool is_sparse() const {
        return false;
    }

    double is_sparse_proportion() const {
        return 0;
    }

    using tatami::Matrix<Value_, Index_>::dense;

    using tatami::Matrix<Value_, Index_>::sparse;

    /********************
     *** Myopic dense ***
     ********************/
private:
    template<bool oracle_, template<bool, bool, typename, typename, class> class Extractor_, typename ... Args_>
    std::unique_ptr<tatami::DenseExtractor<oracle_, Value_, Index_> > raw_dense_internal(bool row, Index_ secondary_length, Args_&& ... args) const {
        if (row) {
            // Remember, the num_chunks_per_column is the number of slabs needed to divide up all the *rows* of the matrix.
            SlabCacheStats stats(my_coordinator.get_chunk_nrow(), secondary_length, my_coordinator.get_num_chunks_per_column(), my_cache_size_in_elements, my_require_minimum_cache);
            if (stats.max_slabs_in_cache > 0) {
                return std::make_unique<Extractor_<false, oracle_, Value_, Index_, Chunk_> >(my_coordinator, stats, row, std::forward<Args_>(args)...);
            } else {
                return std::make_unique<Extractor_<true, oracle_, Value_, Index_, Chunk_> >(my_coordinator, stats, row, std::forward<Args_>(args)...);
            }
        } else {
            // Remember, the num_chunks_per_row is the number of slabs needed to divide up all the *columns* of the matrix.
            SlabCacheStats stats(my_coordinator.get_chunk_ncol(), secondary_length, my_coordinator.get_num_chunks_per_row(), my_cache_size_in_elements, my_require_minimum_cache);
            if (stats.max_slabs_in_cache > 0) {
                return std::make_unique<Extractor_<false, oracle_, Value_, Index_, Chunk_> >(my_coordinator, stats, row, std::forward<Args_>(args)...);
            } else {
                return std::make_unique<Extractor_<true, oracle_, Value_, Index_, Chunk_> >(my_coordinator, stats, row, std::forward<Args_>(args)...);
            }
        }
    }

    template<bool oracle_>
    std::unique_ptr<tatami::DenseExtractor<oracle_, Value_, Index_> > dense_internal(bool row, tatami::MaybeOracle<oracle_, Index_> oracle, const tatami::Options&) const {
        auto secondary = (row ? my_coordinator.get_ncol() : my_coordinator.get_nrow());
        return raw_dense_internal<oracle_, CustomChunkedMatrix_internal::DenseFull>(row, secondary, std::move(oracle));
    }

    template<bool oracle_>
    std::unique_ptr<tatami::DenseExtractor<oracle_, Value_, Index_> > dense_internal(
        bool row, 
        tatami::MaybeOracle<oracle_, Index_> oracle, 
        Index_ block_start, 
        Index_ block_length, 
        const tatami::Options&) 
    const {
        return raw_dense_internal<oracle_, CustomChunkedMatrix_internal::DenseBlock>(row, block_length, std::move(oracle), block_start, block_length);
    }

    template<bool oracle_>
    std::unique_ptr<tatami::DenseExtractor<oracle_, Value_, Index_> > dense_internal(
        bool row, 
        tatami::MaybeOracle<oracle_, Index_> oracle, 
        tatami::VectorPtr<Index_> indices_ptr, 
        const tatami::Options&) 
    const {
        auto num_indices = indices_ptr->size();
        return raw_dense_internal<oracle_, CustomChunkedMatrix_internal::DenseIndex>(row, num_indices, std::move(oracle), std::move(indices_ptr));
    }

public:
    std::unique_ptr<tatami::MyopicDenseExtractor<Value_, Index_> > dense(bool row, const tatami::Options& opt) const {
        return dense_internal<false>(row, false, opt);
    }

    std::unique_ptr<tatami::MyopicDenseExtractor<Value_, Index_> > dense(bool row, Index_ block_start, Index_ block_length, const tatami::Options& opt) const {
        return dense_internal<false>(row, false, block_start, block_length, opt);
    }

    std::unique_ptr<tatami::MyopicDenseExtractor<Value_, Index_> > dense(bool row, tatami::VectorPtr<Index_> indices_ptr, const tatami::Options& opt) const {
        return dense_internal<false>(row, false, std::move(indices_ptr), opt);
    }

    /**********************
     *** Oracular dense ***
     **********************/
public:
    std::unique_ptr<tatami::OracularDenseExtractor<Value_, Index_> > dense(
        bool row, 
        std::shared_ptr<const tatami::Oracle<Index_> > oracle, 
        const tatami::Options& opt) 
    const {
        return dense_internal<true>(row, std::move(oracle), opt);
    }

    std::unique_ptr<tatami::OracularDenseExtractor<Value_, Index_> > dense(
        bool row, 
        std::shared_ptr<const tatami::Oracle<Index_> > oracle, 
        Index_ block_start, 
        Index_ block_length, 
        const tatami::Options& opt) 
    const {
        return dense_internal<true>(row, std::move(oracle), block_start, block_length, opt);
    }

    std::unique_ptr<tatami::OracularDenseExtractor<Value_, Index_> > dense(
        bool row, 
        std::shared_ptr<const tatami::Oracle<Index_> > oracle, 
        tatami::VectorPtr<Index_> indices_ptr, 
        const tatami::Options& opt) 
    const {
        return dense_internal<true>(row, std::move(oracle), std::move(indices_ptr), opt);
    }

    /*********************
     *** Myopic sparse ***
     *********************/
public:
    std::unique_ptr<tatami::MyopicSparseExtractor<Value_, Index_> > sparse(bool row, const tatami::Options& opt) const {
        return std::make_unique<tatami::FullSparsifiedWrapper<false, Value_, Index_> >(dense(row, opt), my_coordinator.get_secondary_dim(row), opt);
    }

    std::unique_ptr<tatami::MyopicSparseExtractor<Value_, Index_> > sparse(bool row, Index_ block_start, Index_ block_length, const tatami::Options& opt) const {
        return std::make_unique<tatami::BlockSparsifiedWrapper<false, Value_, Index_> >(dense(row, block_start, block_length, opt), block_start, block_length, opt);
    }

    std::unique_ptr<tatami::MyopicSparseExtractor<Value_, Index_> > sparse(bool row, tatami::VectorPtr<Index_> indices_ptr, const tatami::Options& opt) const {
        auto d = dense(row, indices_ptr, opt);
        return std::make_unique<tatami::IndexSparsifiedWrapper<false, Value_, Index_> >(std::move(d), std::move(indices_ptr), opt);
    }

    /***********************
     *** Oracular sparse ***
     ***********************/
public:
    std::unique_ptr<tatami::OracularSparseExtractor<Value_, Index_> > sparse(
        bool row, 
        std::shared_ptr<const tatami::Oracle<Index_> > oracle, 
        const tatami::Options& opt)
    const {
        return std::make_unique<tatami::FullSparsifiedWrapper<true, Value_, Index_> >(dense(row, std::move(oracle), opt), my_coordinator.get_secondary_dim(row), opt);
    }

    std::unique_ptr<tatami::OracularSparseExtractor<Value_, Index_> > sparse(
        bool row, 
        std::shared_ptr<const tatami::Oracle<Index_> > oracle, 
        Index_ block_start, 
        Index_ block_length, 
        const tatami::Options& opt) 
    const {
        return std::make_unique<tatami::BlockSparsifiedWrapper<true, Value_, Index_> >(dense(row, std::move(oracle), block_start, block_length, opt), block_start, block_length, opt);
    }

    std::unique_ptr<tatami::OracularSparseExtractor<Value_, Index_> > sparse(
        bool row, 
        std::shared_ptr<const tatami::Oracle<Index_> > oracle, 
        tatami::VectorPtr<Index_> indices_ptr, 
        const tatami::Options& opt)
    const {
        auto d = dense(row, std::move(oracle), indices_ptr, opt);
        return std::make_unique<tatami::IndexSparsifiedWrapper<true, Value_, Index_> >(std::move(d), std::move(indices_ptr), opt);
    }
};

}

#endif
