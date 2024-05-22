#ifndef TATAMI_R_UNKNOWNMATRIX_HPP
#define TATAMI_R_UNKNOWNMATRIX_HPP

#include "Rcpp.h"
#include "tatami/tatami.hpp"
#include "dense_extractor.hpp"
#include "sparse_extractor.hpp"

#include <vector>
#include <memory>
#include <string>
#include <stdexcept>

#include "parallelize.hpp"

namespace tatami_r {

/**
 * @brief Options for data extraction from an `UnknownMatrix`.
 */
struct UnknownMatrixOptions {
    /**
     * Size of the cache, in bytes.
     * If -1, this is determined from `DelayedArray::getAutoBlockSize()`.
     */
    size_t maximum_cache_size = -1;

    /**
     * Whether to automatically enforce a minimum size for the cache, regardless of `maximum_cache_size`.
     * This minimum is chosen to ensure that all chunks overlapping one row (or a slice/subset thereof) can be retained in memory,
     * so that the same chunks are not repeatedly re-read from disk when iterating over consecutive rows/columns of the matrix.
     */
    bool require_minimum_cache = true;
};

/**
 * @brief Unknown matrix-like object in R.
 *
 * @tparam Value_ Numeric type of data value for the interface.
 * @tparam Index_ Integer type for the row/column indices, for the interface.
 *
 * Pull data out of an unknown matrix-like object by calling methods from the [**DelayedArray**](https://bioconductor.org/packages/DelayedArray) package via **Rcpp**.
 * This effectively extends **tatami** to work with any abstract numeric matrix that might be consumed by an R function.
 * 
 * Instances of class should only be constructed and destroyed in a serial context, specifically on the same thread running R itself. 
 * Calls to its methods may be parallelized but some additional effort is required to serialize calls to the R API; see `executor()` for more details.
 */
template<typename Value_, typename Index_, typename CachedValue_ = Value_, typename CachedIndex_ = Index_>
class UnknownMatrix : public tatami::Matrix<Value_, Index_> {
public:
    /**
     * This constructor should only be called in a serial context, as the (default) construction of **Rcpp** objects may call the R API.
     *
     * @param seed A matrix-like R object.
     * @param opt Extraction options.
     */
    UnknownMatrix(Rcpp::RObject seed, const UnknownMatrixOptions& opt) : 
        my_original_seed(seed), 
        my_delayed_env(Rcpp::Environment::namespace_env("DelayedArray")),
        my_sparse_env(Rcpp::Environment::namespace_env("SparseArray")),
        my_dense_extractor(my_delayed_env["extract_array"]),
        my_sparse_extractor(my_sparse_env["extract_sparse_array"])
    {
        // We assume the constructor only occurs on the main thread, so we
        // won't bother locking things up. I'm also not sure that the
        // operations in the initialization list are thread-safe.

        {
            auto base = Rcpp::Environment::base_env();
            Rcpp::Function fun = base["dim"];
            Rcpp::RObject output = fun(seed);
            if (output.sexp_type() != INTSXP) {
                auto ctype = get_class_name(my_original_seed);
                throw std::runtime_error("'dim(<" + ctype + ">)' should return an integer vector");
            }
            Rcpp::IntegerVector dims(output);
            if (dims.size() != 2 || dims[0] < 0 || dims[1] < 0) {
                auto ctype = get_class_name(my_original_seed);
                throw std::runtime_error("'dim(<" + ctype + ">)' should contain two non-negative integers");
            }

            my_nrow = dims[0];
            my_ncol = dims[1];
        }

        {
            Rcpp::Function fun = my_delayed_env["is_sparse"];
            Rcpp::LogicalVector is_sparse = fun(seed);
            if (is_sparse.size() != 1) {
                auto ctype = get_class_name(my_original_seed);
                throw std::runtime_error("'is_sparse(<" + ctype + ">)' should return a logical vector of length 1");
            }
            my_sparse = (is_sparse[0] != 0);
        }

        {
            my_row_chunk_map.resize(my_nrow);
            my_col_chunk_map.resize(my_ncol);

            Rcpp::Function fun = my_delayed_env["chunkGrid"];
            Rcpp::RObject grid = fun(seed);

            if (grid == R_NilValue) {
                my_row_max_chunk_size = 1;
                my_col_max_chunk_size = 1;
                std::iota(my_row_chunk_map.begin(), my_row_chunk_map.end(), static_cast<Index_>(0));
                std::iota(my_col_chunk_map.begin(), my_col_chunk_map.end(), static_cast<Index_>(0));
                my_row_chunk_ticks.resize(my_nrow + 1);
                std::iota(my_row_chunk_ticks.begin(), my_row_chunk_ticks.end(), static_cast<Index_>(0));
                my_col_chunk_ticks.resize(my_ncol + 1);
                std::iota(my_col_chunk_ticks.begin(), my_col_chunk_ticks.end(), static_cast<Index_>(0));

                // Both dense and sparse inputs are implicitly column-major, so
                // if there isn't chunking information to the contrary, we'll
                // favor extraction of the columns.
                my_prefer_rows = false;

            } else {
                auto grid_cls = get_class_name(grid);

                if (grid_cls == "RegularArrayGrid") {
                    Rcpp::IntegerVector spacings(Rcpp::RObject(grid.slot("spacings")));
                    if (spacings.size() != 2) {
                        auto ctype = get_class_name(seed);
                        throw std::runtime_error("'chunkGrid(<" + ctype + ">)@spacings' should be an integer vector of length 2 with non-negative values");
                    }

                    auto populate = [](Index_ extent, Index_ spacing, std::vector<Index_>& map, std::vector<Index_>& ticks) {
                        if (spacing == 0) {
                            ticks.push_back(0);
                        } else {
                            ticks.reserve((extent / spacing) + (extent % spacing > 0) + 1);
                            Index_ start = 0;
                            ticks.push_back(start);
                            while (start != extent) {
                                auto to_fill = std::min(spacing, extent - start);
                                std::fill_n(map.begin() + start, to_fill, ticks.size() - 1);
                                start += to_fill;
                                ticks.push_back(start);
                            }
                        }
                    };

                    my_row_max_chunk_size = spacings[0];
                    populate(my_nrow, my_row_max_chunk_size, my_row_chunk_map, my_row_chunk_ticks);
                    my_col_max_chunk_size = spacings[1];
                    populate(my_ncol, my_col_max_chunk_size, my_col_chunk_map, my_col_chunk_ticks);

                } else if (grid_cls == "ArbitraryArrayGrid") {
                    Rcpp::List ticks(Rcpp::RObject(grid.slot("tickmarks")));
                    if (ticks.size() != 2) {
                        auto ctype = get_class_name(seed);
                        throw std::runtime_error("'chunkGrid(<" + ctype + ">)@tickmarks' should return a list of length 2");
                    }

                    auto populate = [](Index_ extent, const Rcpp::IntegerVector& ticks, std::vector<Index_>& map, std::vector<Index_>& new_ticks, Index_& max_chunk_size) {
                        if (ticks.size() == 0 || ticks[ticks.size() - 1] != static_cast<int>(extent)) {
                            throw std::runtime_error("invalid ticks returned by 'chunkGrid'");
                        }
                        new_ticks.resize(ticks.size() + 1);
                        std::copy(ticks.begin(), ticks.end(), new_ticks.begin() + 1);

                        max_chunk_size = 0;
                        int start = 0;
                        map.resize(extent);
                        Index_ counter = 0;

                        for (auto t : ticks) {
                            if (t < start) {
                                throw std::runtime_error("invalid ticks returned by 'chunkGrid'");
                            }
                            Index_ to_fill = t - start;
                            if (to_fill > max_chunk_size) {
                                max_chunk_size = to_fill;
                            }
                            std::fill_n(map.begin() + start, to_fill, counter);
                            ++counter;
                            start = t;
                        }
                    };

                    Rcpp::IntegerVector first(ticks[0]);
                    populate(my_nrow, first, my_row_chunk_map, my_row_chunk_ticks, my_row_max_chunk_size);
                    Rcpp::IntegerVector second(ticks[1]);
                    populate(my_ncol, second, my_col_chunk_map, my_col_chunk_ticks, my_col_max_chunk_size);

                } else {
                    auto ctype = get_class_name(seed);
                    throw std::runtime_error("instance of unknown class '" + grid_cls + "' returned by 'chunkGrid(<" + ctype + ">)");
                }

                // Choose the dimension that requires pulling out fewer chunks.
                auto chunks_per_row = my_col_chunk_ticks.size() - 1;
                auto chunks_per_col = my_row_chunk_ticks.size() - 1;
                my_prefer_rows = chunks_per_row <= chunks_per_col;
            }
        }

        my_cache_size_in_bytes = opt.maximum_cache_size;
        my_require_minimum_cache = opt.require_minimum_cache;
        if (my_cache_size_in_bytes == static_cast<size_t>(-1)) {
            Rcpp::Function fun = my_delayed_env["getAutoBlockSize"];
            Rcpp::NumericVector bsize = fun();
            if (bsize.size() != 1 || bsize[0] < 0) {
                throw std::runtime_error("'getAutoBlockSize()' should return a non-negative number of bytes");
            }
            my_cache_size_in_bytes = bsize[0];
        }
    }

    /**
     * This constructor overload uses the default `Options()`.
     * Again, it should only be called in a serial context, as the (default) construction of **Rcpp** objects may call the R API.
     *
     * @param seed A matrix-like R object.
     */
    UnknownMatrix(Rcpp::RObject seed) : UnknownMatrix(std::move(seed), UnknownMatrixOptions()) {}

private:
    Index_ my_nrow, my_ncol;
    bool my_sparse, my_prefer_rows;

    std::vector<Index_> my_row_chunk_map, my_col_chunk_map;
    std::vector<Index_> my_row_chunk_ticks, my_col_chunk_ticks;

    // To decide how many chunks to store in the cache, we pretend the largest
    // chunk is a good representative. This is a bit suboptimal for irregular
    // chunks but the LruSlabCache class doesn't have a good way of dealing
    // with this right now. The fundamental problem is that variable slabs will
    // either (i) all reach the maximum allocation eventually, if slabs are
    // reused, or (ii) require lots of allocations, if slabs are not reused, or
    // (iii) require manual defragmentation, if slabs are reused in a manner
    // that avoids inflation to the maximum allocation.
    Index_ my_row_max_chunk_size, my_col_max_chunk_size;

    size_t my_cache_size_in_bytes;
    bool my_require_minimum_cache;

    Rcpp::RObject my_original_seed;
    Rcpp::Environment my_delayed_env, my_sparse_env;
    Rcpp::Function my_dense_extractor, my_sparse_extractor;

public:
    Index_ nrow() const {
        return my_nrow;
    }

    Index_ ncol() const {
        return my_ncol;
    }

    bool is_sparse() const {
        return my_sparse;
    }

    double is_sparse_proportion() const {
        return static_cast<double>(my_sparse);
    }

    bool prefer_rows() const {
        return my_prefer_rows;
    }

    double prefer_rows_proportion() const {
        return static_cast<double>(my_prefer_rows);
    }

    bool uses_oracle(bool) const {
        return true;
    }

private:
    Index_ max_primary_chunk_length(bool row) const {
        return (row ? my_row_max_chunk_size : my_col_max_chunk_size);
    }

    Index_ secondary_dim(bool row) const {
        return (row ? my_ncol : my_nrow);
    }

    const std::vector<Index_>& chunk_ticks(bool row) const {
        if (row) {
            return my_row_chunk_ticks;
        } else {
            return my_col_chunk_ticks;
        }
    }

    const std::vector<Index_>& chunk_map(bool row) const {
        if (row) {
            return my_row_chunk_map;
        } else {
            return my_col_chunk_map;
        }
    }

    /********************
     *** Myopic dense ***
     ********************/
private:
    template<
        bool oracle_, 
        template <bool, bool, typename, typename, typename> class FromDense_,
        template <bool, bool, typename, typename, typename, typename> class FromSparse_,
        typename ... Args_
    >
    std::unique_ptr<tatami::DenseExtractor<oracle_, Value_, Index_> > populate_dense_internal(bool row, Index_ non_target_length, tatami::MaybeOracle<oracle_, Index_> oracle, Args_&& ... args) const {
        std::unique_ptr<tatami::DenseExtractor<oracle_, Value_, Index_> > output;

        Index_ max_target_chunk_length = max_primary_chunk_length(row);
        tatami_chunked::SlabCacheStats stats(max_target_chunk_length, non_target_length, my_cache_size_in_bytes, sizeof(CachedValue_), my_require_minimum_cache);

        const auto& map = chunk_map(row);
        const auto& ticks = chunk_ticks(row);
        bool solo = (stats.max_slabs_in_cache == 0);

#ifdef TATAMI_R_PARALLELIZE_UNKNOWN 
        // This involves some Rcpp initializations, so we lock it just in case.
        auto& mexec = executor();
        mexec.run([&]() -> void {
#endif

        if (!my_sparse) {
            if (solo) {
                typedef FromDense_<true, oracle_, Value_, Index_, CachedValue_> ShortDense;
                output.reset(new ShortDense(my_original_seed, my_dense_extractor, row, std::move(oracle), std::forward<Args_>(args)..., ticks, map, stats));
            } else {
                typedef FromDense_<false, oracle_, Value_, Index_, CachedValue_> ShortDense;
                output.reset(new ShortDense(my_original_seed, my_dense_extractor, row, std::move(oracle), std::forward<Args_>(args)..., ticks, map, stats));
            }
        } else {
            if (solo) {
                typedef FromSparse_<true, oracle_, Value_, Index_, CachedValue_, CachedIndex_> ShortSparse;
                output.reset(new ShortSparse(my_original_seed, my_sparse_extractor, row, std::move(oracle), std::forward<Args_>(args)..., max_target_chunk_length, ticks, map, stats));
            } else {
                typedef FromSparse_<false, oracle_, Value_, Index_, CachedValue_, CachedIndex_> ShortSparse;
                output.reset(new ShortSparse(my_original_seed, my_sparse_extractor, row, std::move(oracle), std::forward<Args_>(args)..., max_target_chunk_length, ticks, map, stats));
            }
        }

#ifdef TATAMI_R_PARALLELIZE_UNKNOWN 
        });
#endif

        return output;
    }

    template<bool oracle_>
    std::unique_ptr<tatami::DenseExtractor<oracle_, Value_, Index_> > populate_dense(bool row, tatami::MaybeOracle<oracle_, Index_> ora, const tatami::Options&) const {
        Index_ non_target_dim = secondary_dim(row);
        return populate_dense_internal<oracle_, UnknownMatrix_internal::DenseFull, UnknownMatrix_internal::DensifiedSparseFull>(row, non_target_dim, std::move(ora), non_target_dim);
    }

    template<bool oracle_>
    std::unique_ptr<tatami::DenseExtractor<oracle_, Value_, Index_> > populate_dense(bool row, tatami::MaybeOracle<oracle_, Index_> ora, Index_ block_start, Index_ block_length, const tatami::Options&) const {
        return populate_dense_internal<oracle_, UnknownMatrix_internal::DenseBlock, UnknownMatrix_internal::DensifiedSparseBlock>(row, block_length, std::move(ora), block_start, block_length);
    }

    template<bool oracle_>
    std::unique_ptr<tatami::DenseExtractor<oracle_, Value_, Index_> > populate_dense(bool row, tatami::MaybeOracle<oracle_, Index_> ora, tatami::VectorPtr<Index_> indices_ptr, const tatami::Options&) const {
        Index_ nidx = indices_ptr->size();
        return populate_dense_internal<oracle_, UnknownMatrix_internal::DenseIndexed, UnknownMatrix_internal::DensifiedSparseIndexed>(row, nidx, std::move(ora), std::move(indices_ptr));
    }

public:
    std::unique_ptr<tatami::MyopicDenseExtractor<Value_, Index_> > dense(bool row, const tatami::Options& opt) const {
        return populate_dense<false>(row, false, opt); 
    }

    std::unique_ptr<tatami::MyopicDenseExtractor<Value_, Index_> > dense(bool row, Index_ block_start, Index_ block_length, const tatami::Options& opt) const {
        return populate_dense<false>(row, false, block_start, block_length, opt); 
    }

    std::unique_ptr<tatami::MyopicDenseExtractor<Value_, Index_> > dense(bool row, tatami::VectorPtr<Index_> indices_ptr, const tatami::Options& opt) const {
        return populate_dense<false>(row, false, std::move(indices_ptr), opt); 
    }

    /**********************
     *** Oracular dense ***
     **********************/
public:
    std::unique_ptr<tatami::OracularDenseExtractor<Value_, Index_> > dense(bool row, std::shared_ptr<const tatami::Oracle<Index_> > ora, const tatami::Options& opt) const {
        return populate_dense<true>(row, std::move(ora), opt); 
    }

    std::unique_ptr<tatami::OracularDenseExtractor<Value_, Index_> > dense(bool row, std::shared_ptr<const tatami::Oracle<Index_> > ora, Index_ block_start, Index_ block_length, const tatami::Options& opt) const {
        return populate_dense<true>(row, std::move(ora), block_start, block_length, opt); 
    }

    std::unique_ptr<tatami::OracularDenseExtractor<Value_, Index_> > dense(bool row, std::shared_ptr<const tatami::Oracle<Index_> > ora, tatami::VectorPtr<Index_> indices_ptr, const tatami::Options& opt) const {
        return populate_dense<true>(row, std::move(ora), std::move(indices_ptr), opt); 
    }

    /*********************
     *** Myopic sparse ***
     *********************/
public:
    template<
        bool oracle_, 
        template<bool, bool, typename, typename, typename, typename> class FromSparse_,
        typename ... Args_
    >
    std::unique_ptr<tatami::SparseExtractor<oracle_, Value_, Index_> > populate_sparse_internal(
        bool row,
        Index_ non_target_length, 
        tatami::MaybeOracle<oracle_, Index_> oracle, 
        const tatami::Options& opt, 
        Args_&& ... args) 
    const {
        std::unique_ptr<tatami::SparseExtractor<oracle_, Value_, Index_> > output;

        Index_ max_target_chunk_length = max_primary_chunk_length(row);
        tatami_chunked::SlabCacheStats stats(
            max_target_chunk_length,
            non_target_length, 
            my_cache_size_in_bytes, 
            (opt.sparse_extract_index ? sizeof(CachedIndex_) : 0) + (opt.sparse_extract_value ? sizeof(CachedValue_) : 0),
            my_require_minimum_cache
        );

        const auto& map = chunk_map(row);
        const auto& ticks = chunk_ticks(row);
        bool needs_value = opt.sparse_extract_value;
        bool needs_index = opt.sparse_extract_index;
        bool solo = stats.max_slabs_in_cache == 0;

#ifdef TATAMI_R_PARALLELIZE_UNKNOWN 
        // This involves some Rcpp initializations, so we lock it just in case.
        auto& mexec = executor();
        mexec.run([&]() -> void {
#endif

        if (solo) {
            typedef FromSparse_<true, oracle_, Value_, Index_, CachedValue_, CachedIndex_> ShortSparse;
            output.reset(new ShortSparse(my_original_seed, my_sparse_extractor, row, std::move(oracle), std::forward<Args_>(args)..., max_target_chunk_length, ticks, map, stats, needs_value, needs_index));
        } else {
            typedef FromSparse_<false, oracle_, Value_, Index_, CachedValue_, CachedIndex_> ShortSparse;
            output.reset(new ShortSparse(my_original_seed, my_sparse_extractor, row, std::move(oracle), std::forward<Args_>(args)..., max_target_chunk_length, ticks, map, stats, needs_value, needs_index));
        }

#ifdef TATAMI_R_PARALLELIZE_UNKNOWN 
        });
#endif

        return output;
    }

    template<bool oracle_>
    std::unique_ptr<tatami::SparseExtractor<oracle_, Value_, Index_> > populate_sparse(bool row, tatami::MaybeOracle<oracle_, Index_> ora, const tatami::Options& opt) const {
        Index_ non_target_dim = secondary_dim(row);
        return populate_sparse_internal<oracle_, UnknownMatrix_internal::SparseFull>(row, non_target_dim, std::move(ora), opt, non_target_dim); 
    }

    template<bool oracle_>
    std::unique_ptr<tatami::SparseExtractor<oracle_, Value_, Index_> > populate_sparse(bool row, tatami::MaybeOracle<oracle_, Index_> ora, Index_ block_start, Index_ block_length, const tatami::Options& opt) const {
        return populate_sparse_internal<oracle_, UnknownMatrix_internal::SparseBlock>(row, block_length, std::move(ora), opt, block_start, block_length);
    }

    template<bool oracle_>
    std::unique_ptr<tatami::SparseExtractor<oracle_, Value_, Index_> > populate_sparse(bool row, tatami::MaybeOracle<oracle_, Index_> ora, tatami::VectorPtr<Index_> indices_ptr, const tatami::Options& opt) const {
        Index_ nidx = indices_ptr->size();
        return populate_sparse_internal<oracle_, UnknownMatrix_internal::SparseIndexed>(row, nidx, std::move(ora), opt, std::move(indices_ptr));
    }

public:
    std::unique_ptr<tatami::MyopicSparseExtractor<Value_, Index_> > sparse(bool row, const tatami::Options& opt) const {
        if (!my_sparse) {
            return std::make_unique<tatami::FullSparsifiedWrapper<false, Value_, Index_> >(dense(row, opt), secondary_dim(row), opt);
        } else {
            return populate_sparse<false>(row, false, opt); 
        }
    }

    std::unique_ptr<tatami::MyopicSparseExtractor<Value_, Index_> > sparse(bool row, Index_ block_start, Index_ block_length, const tatami::Options& opt) const {
        if (!my_sparse) {
            return std::make_unique<tatami::BlockSparsifiedWrapper<false, Value_, Index_> >(dense(row, block_start, block_length, opt), block_start, block_length, opt);
        } else {
            return populate_sparse<false>(row, false, block_start, block_length, opt); 
        }
    }

    std::unique_ptr<tatami::MyopicSparseExtractor<Value_, Index_> > sparse(bool row, tatami::VectorPtr<Index_> indices_ptr, const tatami::Options& opt) const {
        if (!my_sparse) {
            auto index_copy = indices_ptr;
            return std::make_unique<tatami::IndexSparsifiedWrapper<false, Value_, Index_> >(dense(row, std::move(indices_ptr), opt), std::move(index_copy), opt);
        } else {
            return populate_sparse<false>(row, false, std::move(indices_ptr), opt); 
        }
    }

    /**********************
     *** Oracular sparse ***
     **********************/
public:
    std::unique_ptr<tatami::OracularSparseExtractor<Value_, Index_> > sparse(bool row, std::shared_ptr<const tatami::Oracle<Index_> > ora, const tatami::Options& opt) const {
        if (!my_sparse) {
            return std::make_unique<tatami::FullSparsifiedWrapper<true, Value_, Index_> >(dense(row, std::move(ora), opt), secondary_dim(row), opt);
        } else {
            return populate_sparse<true>(row, std::move(ora), opt); 
        }
    }

    std::unique_ptr<tatami::OracularSparseExtractor<Value_, Index_> > sparse(bool row, std::shared_ptr<const tatami::Oracle<Index_> > ora, Index_ block_start, Index_ block_length, const tatami::Options& opt) const {
        if (!my_sparse) {
            return std::make_unique<tatami::BlockSparsifiedWrapper<true, Value_, Index_> >(dense(row, std::move(ora), block_start, block_length, opt), block_start, block_length, opt);
        } else {
            return populate_sparse<true>(row, std::move(ora), block_start, block_length, opt); 
        }
    }

    std::unique_ptr<tatami::OracularSparseExtractor<Value_, Index_> > sparse(bool row, std::shared_ptr<const tatami::Oracle<Index_> > ora, tatami::VectorPtr<Index_> indices_ptr, const tatami::Options& opt) const {
        if (!my_sparse) {
            auto index_copy = indices_ptr;
            return std::make_unique<tatami::IndexSparsifiedWrapper<true, Value_, Index_> >(dense(row, std::move(ora), std::move(indices_ptr), opt), std::move(index_copy), opt);
        } else {
            return populate_sparse<true>(row, std::move(ora), std::move(indices_ptr), opt); 
        }
    }
};

}

#endif
