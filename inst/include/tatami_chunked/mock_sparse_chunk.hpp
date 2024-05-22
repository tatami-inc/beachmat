#ifndef TATAMI_CHUNKED_MOCK_SPARSE_CHUNK_HPP
#define TATAMI_CHUNKED_MOCK_SPARSE_CHUNK_HPP

#include <vector>
#include <cstdint>

/**
 * @file mock_sparse_chunk.hpp
 * @brief Sparse chunk interface to use in a `CustomSparseChunkedMatrix`.
 */

namespace tatami_chunked {

/**
 * @cond
 */
namespace MockSparseChunk_internal {

template<typename InflatedValue_, typename InflatedIndex_>
struct Workspace {
    // Standard compressed sparse members:
    std::vector<InflatedValue_> values;
    std::vector<InflatedIndex_> indices;
    std::vector<size_t> indptrs;

    // Allocation to allow for O(1) mapping of requested indices to sparse indices.
    // This mimics what is done in the indexed sparse extractors in tatami proper.
    std::vector<uint8_t> remap;
};

template<class Blob_>
class Core {
public:
    Core() = default;

    Core(Blob_ chunk) : my_chunk(std::move(chunk)) {}

private:
    Blob_ my_chunk;

    typedef typename Blob_::value_type value_type;

    typedef typename Blob_::index_type index_type;

public:
    auto get_target_chunkdim(bool row) const {
        if (row) {
            return my_chunk.nrow();
        } else {
            return my_chunk.ncol();
        }
    }

    auto get_non_target_chunkdim(bool row) const {
        if (row) {
            return my_chunk.ncol();
        } else {
            return my_chunk.nrow();
        }
    }

    typedef Workspace<value_type, index_type> MyWorkspace;

private:
    template<typename Index_>
    static void refine_start_and_end(size_t& start, size_t& end, Index_ desired_start, Index_ desired_end, Index_ max_end, const std::vector<index_type>& indices) {
        if (desired_start) {
            auto it = indices.begin();
            // Using custom comparator to ensure that we cast to Index_ for signedness-safe comparisons.
            start = std::lower_bound(it + start, it + end, desired_start, [](Index_ a, Index_ b) -> bool { return a < b; }) - it;
        }

        if (desired_end != max_end) {
            if (desired_end == desired_start + 1) {
                if (start != end && static_cast<Index_>(indices[start]) == desired_start) {
                    end = start + 1;
                } else {
                    end = start;
                }
            } else {
                auto it = indices.begin();
                end = std::lower_bound(it + start, it + end, desired_end, [](Index_ a, Index_ b) -> bool { return a < b; }) - it;
            }
        }
    }

    // Building a present/absent mapping for the requested indices to allow O(1) look-up from the sparse matrix indices.
    template<typename Index_>
    static void configure_remap(MyWorkspace& work, const std::vector<Index_>& indices, size_t full) {
        work.remap.resize(full);
        for (auto i : indices) {
            work.remap[i] = 1;
        }
    }

    // Resetting just the affected indices so we can avoid a fill operation over the entire array.
    template<typename Index_>
    static void reset_remap(MyWorkspace& work, const std::vector<Index_>& indices) {
        for (auto i : indices) {
            work.remap[i] = 0;
        }
    }

private:
    template<bool is_block_, typename Index_>
    void fill_target(
        Index_ p, 
        Index_ non_target_start, 
        Index_ non_target_end, 
        Index_ non_target_chunkdim,
        MyWorkspace& work, 
        const std::vector<value_type*>& output_values, 
        const std::vector<Index_*>& output_indices,
        Index_* output_number,
        Index_ shift)
    const {
        size_t start = work.indptrs[p], end = work.indptrs[p + 1];
        if (start >= end) {
            return;
        }

        refine_start_and_end(start, end, non_target_start, non_target_end, non_target_chunkdim, work.indices);

        auto& current_number = output_number[p];
        const bool needs_value = !output_values.empty();
        auto vptr = needs_value ? output_values[p] + current_number : NULL;
        const bool needs_index = !output_indices.empty();
        auto iptr = needs_index ? output_indices[p] + current_number : NULL;

        if constexpr(is_block_) {
            if (needs_value) {
                std::copy(work.values.begin() + start, work.values.begin() + end, vptr);
            }
            if (needs_index) {
                for (size_t i = start; i < end; ++i, ++iptr) {
                    *iptr = work.indices[i] + shift;
                }
            }
            current_number += end - start;

        } else {
            // Assumes that work.remap has been properly configured, see configure_remap().
            for (size_t i = start; i < end; ++i) {
                Index_ target = work.indices[i];
                if (work.remap[target]) {
                    if (needs_value) {
                        *vptr = work.values[i];
                        ++vptr;
                    }
                    if (needs_index) {
                        *iptr = target + shift;
                        ++iptr;
                    }
                    ++current_number;
                }
            }
        }
    }

    template<bool is_block_, typename Index_>
    void fill_secondary(
        Index_ s,
        Index_ target_start, 
        Index_ target_end, 
        Index_ target_chunkdim,
        MyWorkspace& work, 
        const std::vector<value_type*>& output_values, 
        const std::vector<Index_*>& output_indices,
        Index_* output_number,
        Index_ shift)
    const {
        auto start = work.indptrs[s], end = work.indptrs[s + 1];
        if (start >= end) {
            return;
        }

        refine_start_and_end(start, end, target_start, target_end, target_chunkdim, work.indices);

        bool needs_value = !output_values.empty();
        bool needs_index = !output_indices.empty();

        if constexpr(is_block_) {
            for (size_t i = start; i < end; ++i) {
                auto p = work.indices[i];
                auto& num = output_number[p];
                if (needs_value) {
                    output_values[p][num] = work.values[i];
                }
                if (needs_index) {
                    output_indices[p][num] = s + shift;
                }
                ++num;
            }

        } else {
            // Assumes that work.remap has been properly configured, see configure_remap().
            for (size_t i = start; i < end; ++i) {
                Index_ target = work.indices[i];
                if (work.remap[target]) {
                    auto& num = output_number[target];
                    if (needs_value) {
                        output_values[target][num] = work.values[i];
                    }
                    if (needs_index) {
                        output_indices[target][num] = s + shift;
                    }
                    ++num;
                }
            }
        }
    }

public:
    template<typename Index_>
    void extract(
        bool row,
        Index_ target_start,
        Index_ target_length,
        Index_ non_target_start,
        Index_ non_target_length,
        MyWorkspace& work,
        const std::vector<value_type*>& output_values, 
        const std::vector<Index_*>& output_indices,
        Index_* output_number,
        Index_ shift)
    const {
        my_chunk.inflate(work.values, work.indices, work.indptrs);
        Index_ target_end = target_start + target_length;
        Index_ non_target_end = non_target_start + non_target_length;

        if (my_chunk.is_csr() == row) {
            Index_ non_target_chunkdim = get_non_target_chunkdim(row);
            for (Index_ p = target_start; p < target_end; ++p) {
                fill_target<true>(p, non_target_start, non_target_end, non_target_chunkdim, work, output_values, output_indices, output_number, shift);
            }
        } else {
            Index_ target_chunkdim = get_target_chunkdim(row);
            for (Index_ s = non_target_start; s < non_target_end; ++s) {
                fill_secondary<true>(s, target_start, target_end, target_chunkdim, work, output_values, output_indices, output_number, shift);
            }
        }
    }

    template<typename Index_>
    void extract(
        bool row,
        Index_ target_start,
        Index_ target_length,
        const std::vector<Index_>& non_target_indices,
        MyWorkspace& work,
        const std::vector<value_type*>& output_values, 
        const std::vector<Index_*>& output_indices,
        Index_* output_number,
        Index_ shift)
    const {
        my_chunk.inflate(work.values, work.indices, work.indptrs);
        Index_ target_end = target_start + target_length;

        if (my_chunk.is_csr() == row) {
            // non_target_indices is guaranteed to be non-empty, see contracts below.
            auto non_target_start = non_target_indices.front();
            auto non_target_end = non_target_indices.back() + 1; // need 1 past end.
            auto non_target_chunkdim = get_non_target_chunkdim(row);

            configure_remap(work, non_target_indices, non_target_chunkdim);
            for (Index_ p = target_start; p < target_end; ++p) {
                fill_target<false>(p, non_target_start, non_target_end, non_target_chunkdim, work, output_values, output_indices, output_number, shift);
            }
            reset_remap(work, non_target_indices);

        } else {
            Index_ target_chunkdim = get_target_chunkdim(row);
            for (auto s : non_target_indices) {
                fill_secondary<true>(s, target_start, target_end, target_chunkdim, work, output_values, output_indices, output_number, shift);
            }
        }
    }

public:
    template<typename Index_>
    void extract(
        bool row,
        const std::vector<Index_>& target_indices,
        Index_ non_target_start,
        Index_ non_target_length,
        MyWorkspace& work,
        const std::vector<value_type*>& output_values, 
        const std::vector<Index_*>& output_indices,
        Index_* output_number,
        Index_ shift)
    const {
        my_chunk.inflate(work.values, work.indices, work.indptrs);
        Index_ non_target_end = non_target_start + non_target_length;

        if (my_chunk.is_csr() == row) {
            Index_ non_target_chunkdim = get_non_target_chunkdim(row);
            for (auto p : target_indices) {
                fill_target<true>(p, non_target_start, non_target_end, non_target_chunkdim, work, output_values, output_indices, output_number, shift);
            }

        } else {
            // target_indices is guaranteed to be non-empty, see contracts below.
            auto target_start = target_indices.front();
            auto target_end = target_indices.back() + 1; // need 1 past end.
            auto target_chunkdim = get_target_chunkdim(row);

            configure_remap(work, target_indices, target_chunkdim);
            for (Index_ s = non_target_start; s < non_target_end; ++s) {
                fill_secondary<false>(s, target_start, target_end, target_chunkdim, work, output_values, output_indices, output_number, shift);
            }
            reset_remap(work, target_indices);
        }
    }

    template<typename Index_>
    void extract(
        bool row,
        const std::vector<Index_>& target_indices,
        const std::vector<Index_>& non_target_indices,
        MyWorkspace& work,
        const std::vector<value_type*>& output_values, 
        const std::vector<Index_*>& output_indices,
        Index_* output_number,
        Index_ shift)
    const {
        my_chunk.inflate(work.values, work.indices, work.indptrs);

        if (my_chunk.is_csr() == row) {
            // non_target_indices is guaranteed to be non-empty, see contracts below.
            auto non_target_start = non_target_indices.front();
            auto non_target_end = non_target_indices.back() + 1; // need 1 past end.
            auto non_target_chunkdim = get_non_target_chunkdim(row);

            configure_remap(work, non_target_indices, non_target_chunkdim);
            for (auto p : target_indices) {
                fill_target<false>(p, non_target_start, non_target_end, non_target_chunkdim, work, output_values, output_indices, output_number, shift);
            }
            reset_remap(work, non_target_indices);

        } else {
            // target_indices is guaranteed to be non-empty, see contracts below.
            auto target_start = target_indices.front();
            auto target_end = target_indices.back() + 1; // need 1 past end.
            auto target_chunkdim = get_target_chunkdim(row);

            configure_remap(work, target_indices, target_chunkdim);
            for (auto s : non_target_indices) {
                fill_secondary<false>(s, target_start, target_end, target_chunkdim, work, output_values, output_indices, output_number, shift);
            }
            reset_remap(work, target_indices);
        }
    }
};

class MockBlob {
public:
    typedef int index_type;
    typedef double value_type;

    bool is_csr() const {
        return true;
    }

private:
    int my_nrow, my_ncol;
    std::vector<double> my_values;
    std::vector<int> my_indices;
    std::vector<size_t> my_pointers;

public:
    MockBlob() = default;

    MockBlob(int nrow, int ncol, std::vector<double> values, std::vector<int> indices, std::vector<size_t> pointers) : 
        my_nrow(nrow), my_ncol(ncol), my_values(std::move(values)), my_indices(std::move(indices)), my_pointers(std::move(pointers)) {}

    int nrow() const {
        return my_nrow;
    }

    int ncol() const {
        return my_ncol;
    }

    void inflate(std::vector<double>& vbuffer, std::vector<int>& ibuffer, std::vector<size_t>& pbuffer) const {
        vbuffer.resize(my_values.size());
        std::copy(my_values.begin(), my_values.end(), vbuffer.begin());
        ibuffer.resize(my_indices.size());
        std::copy(my_indices.begin(), my_indices.end(), ibuffer.begin());
        pbuffer.resize(my_pointers.size());
        std::copy(my_pointers.begin(), my_pointers.end(), pbuffer.begin());
    }
};

}
/**
 * @endcond
 */

/**
 * @brief Mock a simple sparse chunk for a `CustomSparseChunkedMatrix`.
 *
 * Mock a simple sparse chunk for use inside a `CustomSparseChunkedMatrix`.
 * Each chunk should represent a 2-dimensional array of numeric values.
 * The interface is "simple" as extraction of any data involves realization of the entire chunk's contents along the target dimension,
 * with no attempt at optimization if only a subset of dimension elements are of interest.
 */
class MockSimpleSparseChunk {
public:
    /**
     * Type of the value stored in this chunk.
     * Implementations can use any numeric type. 
     */
    typedef double value_type;

    /**
     * Temporary workspace for extracting data from the chunk.
     * One instance of this workspace will be re-used in multiple `extract()` calls for the same or even different chunks.
     * Implementations may use any data structure here.
     */
    struct Workspace {
        /**
         * @cond
         */
        // Hiding this here to avoid giving the impression that we NEED to implement this.
        MockSparseChunk_internal::Workspace<value_type, int> work;
        /**
         * @endcond
         */
    };

    /**
     * Whether to extract a subset of elements on the target dimension.
     * This should be set to `false`, otherwise a `MockSubsettedSparseChunk` is expected.
     */
    static constexpr bool use_subset = false;

public:
    /**
     * @cond
     */
    // You can construct this however you like, I don't care.
    MockSimpleSparseChunk() = default;
    MockSimpleSparseChunk(int nr, int nc, std::vector<double> x, std::vector<int> i, std::vector<size_t> p) : 
        core(MockSparseChunk_internal::MockBlob(nr, nc, std::move(x), std::move(i), std::move(p))) {}
    /**
     * @endcond
     */

private:
    MockSparseChunk_internal::Core<MockSparseChunk_internal::MockBlob> core;

public:
    /**
     * Extract all elements of the target dimension into output buffers. 
     * For each element, this method will extract a contiguous block of the non-target dimension.
     *
     * @tparam Index_ Integer type for the row/column indices of the `CustomSparseChunkedMatrix`.
     *
     * @param row Whether to extract rows from the chunk, i.e., the rows are the target dimension.
     * @param non_target_start Index of the start of a contiguous block of the non-target dimension to be extracted.
     * If `row = true`, this is the first column, otherwise it is the first row.
     * @param non_target_length Length of the contiguous block of the non-target dimension to be extracted.
     * If `row = true`, this is the number of columns, otherwise it is the number of rows.
     * This is guaranteed to be positive.
     * @param work Re-usable workspace for extraction from one or more chunks.
     * @param[out] output_values Vector of pointers in which to store the values of non-zero elements.
     * This has length equal to the extent of the target dimension for this chunk.
     * Each pointer corresponds to an element of the target dimension and refers to an array of length no less than the extent of the non-target dimension of the chunk.
     * Alternatively, this vector may be empty, in which case no values should be stored.
     * @param[out] output_indices Vector of vectors in which to store the indices of the non-zero elements along the non-target dimension.
     * This has length equal to the extent of the target dimension for this chunk.
     * Each pointer corresponds to an element of the target dimension and refers to an array of length no less than the extent of the non-target dimension of the chunk.
     * Alternatively, this vector may be empty, in which case no indices should be stored.
     * @param[in,out] output_number Pointer to an array of length equal to the extent of the target dimension.
     * Each entry `i` specifies the number of non-zero elements that are already present in `output_values[i]` and `output_indices[i]`.
     * @param shift Shift to be added to the chunk's reported indices when storing them in `output_indices`.
     *
     * If `row = true`, we would extract data for all rows and a block of columns `[non_target_start, non_target_start + non_target_length)`;
     * conversely, if `row = false`, we would extract data for all columns and a block of rows.
     * Given a target dimension element `p`, the values of non-zero elements from the requested non-target block should be stored at `output_values[p] + output_number[p]`.
     * The non-target indices for those non-zero elements should be increased by `shift` and stored at `output_indices[p] + output_number[p]` in ascending order.
     * `output_number[p]` should then be increased by the number of non-zero entries stored in this manner.
     * This layout allows concatenation of multiple sparse chunks into a single set of vectors for easier fetching in the `CustomSparseChunkedMatrix`.
     *
     * Note that implementions of this method do not need to have the exact same template arguments as shown here, as long as the types can be deduced.
     */
    template<typename Index_>
    void extract(
        bool row,
        Index_ non_target_start, 
        Index_ non_target_length, 
        Workspace& work, 
        const std::vector<value_type*>& output_values,
        const std::vector<Index_*>& output_indices,
        Index_* output_number,
        Index_ shift)
    const {
        core.template extract<Index_>(
            row, 
            0, 
            core.get_target_chunkdim(row), 
            non_target_start, 
            non_target_length, 
            work.work, 
            output_values, 
            output_indices, 
            output_number,
            shift
        );
    }

    /**
     * Extract all elements of the target dimension into output buffers.
     * For each element, this method will extract an indexed subset of the non-target dimension.
     *
     * @tparam Index_ Integer type for the row/column indices of the `CustomSparseChunkedMatrix`.
     *
     * @param row Whether to extract rows from the chunk, i.e., the rows are the target dimension.
     * @param non_target_indices Indices of the elements on the non-target dimension to be extracted.
     * If `row = true`, these are column indices, otherwise these are row indices.
     * This is guaranteed to be non-empty with unique and sorted indices.
     * @param work Re-usable workspace for extraction from one or more chunks.
     * @param[out] output_values Vector of pointers in which to store the values of non-zero elements.
     * This has length equal to the extent of the target dimension for this chunk.
     * Each pointer corresponds to an element of the target dimension and refers to an array of length no less than the extent of the non-target dimension of the chunk.
     * Alternatively, this vector may be empty, in which case no values should be stored.
     * @param[out] output_indices Vector of vectors in which to store the indices of the non-zero elements along the non-target dimension.
     * This has length equal to the extent of the target dimension for this chunk.
     * Each pointer corresponds to an element of the target dimension and refers to an array of length no less than the extent of the non-target dimension of the chunk.
     * Alternatively, this vector may be empty, in which case no indices should be stored.
     * @param[in,out] output_number Pointer to an array of length equal to the extent of the target dimension.
     * Each entry `i` specifies the number of non-zero elements that are already present in `output_values[i]` and `output_indices[i]`.
     * @param shift Shift to be added to the chunk's reported indices when storing them in `output_indices`.
     *
     * If `row = true`, we would extract all rows and a subset of columns in `non_target_indices`;
     * conversely, if `row = false`, we would extract data for all columns and a subset of rows.
     * Given a target dimension element `p`, the values of non-zero elements from the requested non-target subset should be stored at `output_values[p] + output_number[p]`.
     * The non-target indices for those non-zero elements should be increased by `shift` and stored at `output_indices[p] + output_number[p]` in ascending order.
     * `output_number[p]` should then be increased by the number of non-zero entries stored in this manner.
     * This layout allows concatenation of multiple sparse chunks into a single set of vectors for easier fetching in the `CustomSparseChunkedMatrix`.
     *
     * Note that implementions of this method do not need to have the exact same template arguments as shown here, as long as the types can be deduced.
     */
    template<typename Index_>
    void extract(
        bool row,
        const std::vector<Index_>& non_target_indices,
        Workspace& work, 
        const std::vector<value_type*>& output_values,
        const std::vector<Index_*>& output_indices,
        Index_* output_number,
        Index_ shift)
    const {
        core.template extract<Index_>(
            row, 
            0, 
            core.get_target_chunkdim(row), 
            non_target_indices, 
            work.work, 
            output_values, 
            output_indices,
            output_number,
            shift
        );
    }
};

/**
 * @brief Create a sparse chunk for a `CustomSparseChunkedMatrix`.
 *
 * Wraps a sparse blob in a simple chunk interface for use inside a `CustomSparseChunkedMatrix`.
 * Each blob should hold a (possibly compressed) 2-dimensional sparse array of numeric values.
 * The wrapper satisfies the `MockSimpleSparseChunk` interface, but is even simpler;
 * extraction of any data involves realization of the entire blob, with no optimization for subsets of interest along either dimension.
 *
 * The `Blob_` class should provide the following:
 *
 * - `typedef value_type`, specifying the type of the value in the array.
 * - `typedef index_type`, specifying the type of the index in the submatrix.
 * - A `nrow() const` method, returning the number of rows in the array as an integer.
 * - A `ncol() const` method, returning the numbr of columns in the array as an integer.
 * - A `bool is_csr() const` method, indicating whether the realized data is in compressed sparse row (CSR) format. 
 * - A `void inflate(std::vector<value_type>& values, std::vector<Index_>& indices, std::vector<size_t>& pointers) const` method that fills `values`, `indices` and `pointers`,
 *   each with the corresponding field for compressed sparse matrix format.
 *   (The type of `Index_` is the same as that used in the `CustomSparseChunkedMatrix`.)
 *   This should constitute a CSR matrix if `is_csr()` returns true, and a CSC matrix otherwise.
 *
 * @tparam Blob_ Class to represent a simple chunk.
 */
template<class Blob_>
class SimpleSparseChunkWrapper {
    /**
     * @cond
     */
public:
    typedef typename Blob_::value_type value_type;

    typedef MockSparseChunk_internal::Workspace<value_type, typename Blob_::index_type> Workspace;

    static constexpr bool use_subset = false;

    SimpleSparseChunkWrapper() = default;

    SimpleSparseChunkWrapper(Blob_ core) : my_core(std::move(core)) {}

private:
    MockSparseChunk_internal::Core<Blob_> my_core;

public:
    template<typename Index_>
    void extract(
        bool row,
        Index_ non_target_start, 
        Index_ non_target_length, 
        Workspace& work, 
        const std::vector<value_type*>& output_values,
        const std::vector<Index_*>& output_indices,
        Index_* output_number,
        Index_ shift)
    const {
        my_core.template extract<Index_>(
            row, 
            0, 
            my_core.get_target_chunkdim(row), 
            non_target_start, 
            non_target_length, 
            work, 
            output_values, 
            output_indices, 
            output_number,
            shift
        );
    }

    template<typename Index_>
    void extract(
        bool row,
        const std::vector<Index_>& non_target_indices,
        Workspace& work, 
        const std::vector<value_type*>& output_values,
        const std::vector<Index_*>& output_indices,
        Index_* output_number,
        Index_ shift)
    const {
        my_core.template extract<Index_>(
            row, 
            0, 
            my_core.get_target_chunkdim(row), 
            non_target_indices, 
            work, 
            output_values, 
            output_indices,
            output_number,
            shift
        );
    }
    /**
     * @endcond
     */
};

/**
 * @brief Mock a subsettable sparse chunk for a `CustomSparseChunkedMatrix`.
 *
 * Mock a subsettable sparse chunk for use inside a `CustomSparseChunkedMatrix`.
 * Each chunk should represent a (possible compressed) 2-dimensional array of numeric values.
 * The interface is smarter as it only extracts elements of interest along the target dimension.
 * The elements of interest may be either a contiguous block or a indexed subset,
 * as predicted for each chunk from the `OracularSubsettedSlabCache`.
 * This provides some opportunities for optimization if the chunk supports partial reads.
 */
class MockSubsettedSparseChunk {
public:
    /**
     * Type of the value stored in this chunk.
     * Implementations can use any numeric type. 
     */
    typedef double value_type;

    /**
     * Workspace for chunk extraction.
     * One instance of this workspace will be re-used in multiple `extract()` calls for the same or even different chunks.
     * Implementations maye use any data structure here.
     */
    struct Workspace {
        /**
         * @cond
         */
        // Hiding this here to avoid giving the impression that we NEED to implement this.
        MockSparseChunk_internal::Workspace<value_type, int> work;
        /**
         * @endcond
         */
    };

    /**
     * Whether to extract a subset of elements on the target dimension.
     * This should be set to `true`, otherwise a `MockSimpleSparseChunk` is expected.
     */
    static constexpr bool use_subset = true;

    /**
     * @cond
     */
    MockSubsettedSparseChunk() = default;
    MockSubsettedSparseChunk(int nrow, int ncol, std::vector<double> values, std::vector<int> indices, std::vector<size_t> pointers) : 
        my_core(MockSparseChunk_internal::MockBlob(nrow, ncol, std::move(values), std::move(indices), std::move(pointers))) {}
    /**
     * @endcond
     */

private:
    MockSparseChunk_internal::Core<MockSparseChunk_internal::MockBlob> my_core;

public:
    /**
     * Extract a contiguous block of the target dimension into an output buffer.
     * For each element, this method will only extract a contiguous block of the non-target dimension.
     *
     * @tparam Index_ Integer type for the row/column indices of the `CustomSparseChunkedMatrix`.
     *
     * @param row Whether to extract rows from the chunk, i.e., the rows are the target dimension.
     * @param target_start Index of the start of a contiguous block of the target dimension to be extracted.
     * If `row = true`, this is the first row, otherwise it is the first column.
     * @param target_length Length of the contiguous block of the target dimension to be extracted.
     * If `row = true`, this is the number of rows, otherwise it is the number of columns.
     * This is guaranteed to be positive.
     * @param non_target_start Index of the start of a contiguous block of the non-target dimension to be extracted.
     * If `row = true`, this is the first column, otherwise it is the first row.
     * @param non_target_length Length of the contiguous block of the non-target dimension to be extracted.
     * If `row = true`, this is the number of columns, otherwise it is the number of rows.
     * This is guaranteed to be positive.
     * @param work Re-usable workspace for extraction from one or more chunks.
     * @param[out] output_values Vector of pointers in which to store the values of non-zero elements.
     * This has length equal to the extent of the target dimension for this chunk.
     * Each pointer corresponds to an element of the target dimension and refers to an array of length no less than the extent of the non-target dimension of the chunk.
     * Alternatively, this vector may be empty, in which case no values should be stored.
     * @param[out] output_indices Vector of vectors in which to store the indices of the non-zero elements along the non-target dimension.
     * This has length equal to the extent of the target dimension for this chunk.
     * Each pointer corresponds to an element of the target dimension and refers to an array of length no less than the extent of the non-target dimension of the chunk.
     * Alternatively, this vector may be empty, in which case no indices should be stored.
     * @param[in,out] output_number Pointer to an array of length equal to the extent of the target dimension.
     * Each entry `i` specifies the number of non-zero elements that are already present in `output_values[i]` and `output_indices[i]`.
     * @param shift Shift to be added to the chunk's reported indices when storing them in `output_indices`.
     *
     * If `row = true`, we would extract a block of rows `[target_start, target_start + length)` and a block of columns `[non_target_start, non_target_start + non_target_length)`;
     * conversely, if `row = false`, we would extract a block of columns from the `target_*` arguments and the block of rows from the `non_target_*`arguments. 
     * Given a target dimension element `p`, the values of non-zero elements from the requested non-target block should be stored at `output_values[p] + output_number[p]`.
     * The non-target indices for those non-zero elements should be increased by `shift` and stored at `output_indices[p] + output_number[p]` in ascending order.
     * `output_number[p]` should then be increased by the number of non-zero entries stored in this manner.
     * This layout allows concatenation of multiple sparse chunks into a single set of vectors for easier fetching in the `CustomSparseChunkedMatrix`.
     *
     * Note that implementions of this method do not need to have the exact same template arguments as shown here, as long as the types can be deduced.
     */
    template<typename Index_>
    void extract(
        bool row,
        Index_ target_start, 
        Index_ target_length, 
        Index_ non_target_start, 
        Index_ non_target_length, 
        Workspace& work, 
        const std::vector<value_type*>& output_values,
        const std::vector<Index_*>& output_indices,
        Index_* output_number,
        Index_ shift)
    const {
        my_core.extract(
            row, 
            target_start, 
            target_length, 
            non_target_start, 
            non_target_length, 
            work.work,
            output_values,
            output_indices,
            output_number,
            shift
        );
    }

    /**
     * Extract a contiguous block of the target dimension into an output buffer.
     * For each element, this method will only extract an indexed subset of the non-target dimension.
     *
     * @tparam Index_ Integer type for the row/column indices of the `CustomSparseChunkedMatrix`.
     *
     * @param row Whether to extract rows from the chunk, i.e., the rows are the target dimension.
     * @param target_start Index of the start of a contiguous block of the target dimension to be extracted.
     * If `row = true`, this is the first row, otherwise it is the first column.
     * @param target_length Length of the contiguous block of the target dimension to be extracted.
     * If `row = true`, this is the number of rows, otherwise it is the number of columns.
     * This is guaranteed to be positive.
     * @param non_target_indices Indices of the elements on the non-target dimension to be extracted.
     * If `row = true`, these are column indices, otherwise these are row indices.
     * This is guaranteed to be non-empty with unique and sorted indices.
     * @param work Re-usable workspace for extraction from one or more chunks.
     * @param[out] output_values Vector of pointers in which to store the values of non-zero elements.
     * This has length equal to the extent of the target dimension for this chunk.
     * Each pointer corresponds to an element of the target dimension and refers to an array of length no less than the extent of the non-target dimension of the chunk.
     * Alternatively, this vector may be empty, in which case no values should be stored.
     * @param[out] output_indices Vector of vectors in which to store the indices of the non-zero elements along the non-target dimension.
     * This has length equal to the extent of the target dimension for this chunk.
     * Each pointer corresponds to an element of the target dimension and refers to an array of length no less than the extent of the non-target dimension of the chunk.
     * Alternatively, this vector may be empty, in which case no indices should be stored.
     * @param[in,out] output_number Pointer to an array of length equal to the extent of the target dimension.
     * Each entry `i` specifies the number of non-zero elements that are already present in `output_values[i]` and `output_indices[i]`.
     * @param shift Shift to be added to the chunk's reported indices when storing them in `output_indices`.
     *
     * If `row = true`, we would extract a block of rows `[target_start, target_start + length)` and a subset of columns in `non_target_indices`;
     * conversely, if `row = false`, we would extract a block of columns from the `target_*` arguments and a subset of rows from `non_target_indices`.
     * Given a target dimension element `p`, the values of non-zero elements from the requested non-target subset should be stored at `output_values[p] + output_number[p]`.
     * The non-target indices for those non-zero elements should be increased by `shift` and stored at `output_indices[p] + output_number[p]` in ascending order.
     * `output_number[p]` should then be increased by the number of non-zero entries stored in this manner.
     * This layout allows concatenation of multiple sparse chunks into a single set of vectors for easier fetching in the `CustomSparseChunkedMatrix`.
     *
     * Note that implementions of this method do not need to have the exact same template arguments as shown here, as long as the types can be deduced.
     */
    template<typename Index_>
    void extract(
        bool row,
        Index_ target_start, 
        Index_ target_length, 
        const std::vector<Index_>& non_target_indices, 
        Workspace& work, 
        const std::vector<value_type*>& output_values,
        const std::vector<Index_*>& output_indices,
        Index_* output_number,
        Index_ shift)
    const {
        my_core.extract(
            row, 
            target_start, 
            target_length, 
            non_target_indices, 
            work.work, 
            output_values, 
            output_indices, 
            output_number,
            shift
        );
    }

public:
    /**
     * Extract an indexed subset of the target dimension into an output buffer.
     * For each element, this method will only extract a contiguous block of the non-target dimension.
     *
     * @tparam Index_ Integer type for the row/column indices of the `CustomSparseChunkedMatrix`.
     *
     * @param row Whether to extract rows from the chunk, i.e., the rows are the target dimension.
     * @param target_indices Indices of the elements on the target dimension to be extracted.
     * If `row = true`, these are row indices, otherwise these are column indices.
     * This is guaranteed to be non-empty with unique and sorted indices.
     * @param non_target_start Index of the start of a contiguous block of the non-target dimension to be extracted.
     * If `row = true`, this is the first column, otherwise it is the first row.
     * @param non_target_length Length of the contiguous block of the non-target dimension to be extracted.
     * If `row = true`, this is the number of columns, otherwise it is the number of rows.
     * This is guaranteed to be positive.
     * @param work Re-usable workspace for extraction from one or more chunks.
     * @param[out] output_values Vector of pointers in which to store the values of non-zero elements.
     * This has length equal to the extent of the target dimension for this chunk.
     * Each pointer corresponds to an element of the target dimension and refers to an array of length no less than the extent of the non-target dimension of the chunk.
     * Alternatively, this vector may be empty, in which case no values should be stored.
     * @param[out] output_indices Vector of vectors in which to store the indices of the non-zero elements along the non-target dimension.
     * This has length equal to the extent of the target dimension for this chunk.
     * Each pointer corresponds to an element of the target dimension and refers to an array of length no less than the extent of the non-target dimension of the chunk.
     * Alternatively, this vector may be empty, in which case no indices should be stored.
     * @param[in,out] output_number Pointer to an array of length equal to the extent of the target dimension.
     * Each entry `i` specifies the number of non-zero elements that are already present in `output_values[i]` and `output_indices[i]`.
     * @param shift Shift to be added to the chunk's reported indices when storing them in `output_indices`.
     *
     * If `row = true`, we would extract a subset of rows in `target_indices` and a block of columns `[non_target_start, non_target_start + non_target_length)`,
     * conversely, if `row = false`, we would extract a subset of columns from the `target_*` arguments and the block of rows from the `non_target_*` arguments.
     * Given a target dimension element `p`, the values of non-zero elements from the requested non-target block should be stored at `output_values[p] + output_number[p]`.
     * The non-target indices for those non-zero elements should be increased by `shift` and stored at `output_indices[p] + output_number[p]` in ascending order.
     * `output_number[p]` should then be increased by the number of non-zero entries stored in this manner.
     * This layout allows concatenation of multiple sparse chunks into a single set of vectors for easier fetching in the `CustomSparseChunkedMatrix`.
     *
     * Note that implementions of this method do not need to have the exact same template arguments as shown here, as long as the types can be deduced.
     */
    template<typename Index_>
    void extract(
        bool row,
        const std::vector<Index_>& target_indices, 
        Index_ non_target_start, 
        Index_ non_target_length, 
        Workspace& work, 
        const std::vector<value_type*>& output_values,
        const std::vector<Index_*>& output_indices,
        Index_* output_number,
        Index_ shift)
    const {
        my_core.extract(
            row, 
            target_indices, 
            non_target_start, 
            non_target_length, 
            work.work, 
            output_values,
            output_indices,
            output_number,
            shift
        );
    }

    /**
     * Extract an indexed subset of the target dimension into an output buffer.
     * For each element, this method will only extract an indexed subset of the non-target dimension.
     *
     * @tparam Index_ Integer type for the row/column indices of the `CustomSparseChunkedMatrix`.
     *
     * @param row Whether to extract rows from the chunk, i.e., the rows are the target dimension.
     * @param target_indices Indices of the elements on the target dimension to be extracted.
     * If `row = true`, these are row indices, otherwise these are column indices.
     * This is guaranteed to be non-empty with unique and sorted indices.
     * @param non_target_indices Indices of the elements on the non-target dimension to be extracted.
     * If `row = true`, these are column indices, otherwise these are row indices.
     * This is guaranteed to be non-empty with unique and sorted indices.
     * @param work Re-usable workspace for extraction from one or more chunks.
     * @param[out] output_values Vector of pointers in which to store the values of non-zero elements.
     * This has length equal to the extent of the target dimension for this chunk.
     * Each pointer corresponds to an element of the target dimension and refers to an array of length no less than the extent of the non-target dimension of the chunk.
     * Alternatively, this vector may be empty, in which case no values should be stored.
     * @param[out] output_indices Vector of vectors in which to store the indices of the non-zero elements along the non-target dimension.
     * This has length equal to the extent of the target dimension for this chunk.
     * Each pointer corresponds to an element of the target dimension and refers to an array of length no less than the extent of the non-target dimension of the chunk.
     * Alternatively, this vector may be empty, in which case no indices should be stored.
     * @param[in,out] output_number Pointer to an array of length equal to the extent of the target dimension.
     * Each entry `i` specifies the number of non-zero elements that are already present in `output_values[i]` and `output_indices[i]`.
     * @param shift Shift to be added to the chunk's reported indices when storing them in `output_indices`.
     *
     * If `row = true`, we would extract a subset of rows in `target_indices` and a subset of columns in `non_target_indices`.
     * conversely, if `row = false`, we would extract a subset of columns in `target_indices` and the subset of rows in `non_target_indices`.
     * Given a target dimension element `p`, the values of non-zero elements from the requested non-target subset should be stored at `output_values[p] + output_number[p]`.
     * The non-target indices for those non-zero elements should be increased by `shift` and stored at `output_indices[p] + output_number[p]` in ascending order.
     * `output_number[p]` should then be increased by the number of non-zero entries stored in this manner.
     * This layout allows concatenation of multiple sparse chunks into a single set of vectors for easier fetching in the `CustomSparseChunkedMatrix`.
     *
     * Note that implementions of this method do not need to have the exact same template arguments as shown here, as long as the types can be deduced.
     */
    template<typename Index_>
    void extract(
        bool row,
        const std::vector<Index_>& target_indices, 
        const std::vector<Index_>& non_target_indices, 
        Workspace& work, 
        const std::vector<value_type*>& output_values,
        const std::vector<Index_*>& output_indices,
        Index_* output_number,
        Index_ shift)
    const {
        my_core.extract(
            row, 
            target_indices, 
            non_target_indices, 
            work.work, 
            output_values, 
            output_indices, 
            output_number,
            shift 
        );
    }
};

}

#endif
