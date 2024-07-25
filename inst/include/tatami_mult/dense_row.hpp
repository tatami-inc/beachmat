#ifndef TATAMI_MULT_DENSE_ROW_HPP
#define TATAMI_MULT_DENSE_ROW_HPP

#include <vector>
#include <cstdint>

#include "tatami/tatami.hpp"
#include "utils.hpp"

namespace tatami_mult {

namespace internal {

template<typename Value_, typename Index_, typename Right_, typename Output_>
void dense_row_vector(const tatami::Matrix<Value_, Index_>& matrix, const Right_* rhs, Output_* output, int num_threads) {
    Index_ NR = matrix.nrow();
    Index_ NC = matrix.ncol();

    tatami::parallelize([&](size_t, Index_ start, Index_ length) {
        auto ext = tatami::consecutive_extractor<false>(&matrix, true, start, length);
        std::vector<Value_> buffer(NC);

        for (Index_ r = start, end = start + length; r < end; ++r) {
            auto ptr = ext->fetch(buffer.data());
            output[r] = std::inner_product(ptr, ptr + NC, rhs, static_cast<Output_>(0));
        }

    }, NR, num_threads);
}

template<typename Value_, typename Index_, typename Right_, typename Output_>
void dense_row_vectors(const tatami::Matrix<Value_, Index_>& matrix, const std::vector<Right_*>& rhs, const std::vector<Output_*>& output, int num_threads) {
    Index_ NR = matrix.nrow();
    Index_ NC = matrix.ncol();
    size_t num_rhs = rhs.size();

    tatami::parallelize([&](size_t, Index_ start, Index_ length) {
        auto ext = tatami::consecutive_extractor<false>(&matrix, true, start, length);
        std::vector<Value_> buffer(NC);

        for (Index_ r = start, end = start + length; r < end; ++r) {
            auto ptr = ext->fetch(buffer.data());
            for (size_t j = 0; j < num_rhs; ++j) {
                output[j][r] = std::inner_product(ptr, ptr + NC, rhs[j], static_cast<Output_>(0));
            }
        }

    }, NR, num_threads);
}

template<typename Value_, typename Index_, typename RightValue_, typename RightIndex_, typename Output_>
void dense_row_tatami_dense(const tatami::Matrix<Value_, Index_>& matrix, const tatami::Matrix<RightValue_, RightIndex_>& rhs, Output_* output, size_t row_shift, size_t col_shift, int num_threads) {
    Index_ NR = matrix.nrow();
    Index_ NC = matrix.ncol();
    RightIndex_ rhs_col = rhs.ncol();

    tatami::parallelize([&](size_t, Index_ start, Index_ length) {
        auto ext = tatami::consecutive_extractor<false>(&matrix, true, start, length);
        std::vector<Value_> buffer(NC);
        std::vector<RightValue_> rbuffer(NC);

        for (Index_ r = start, end = start + length; r < end; ++r) {
            auto ptr = ext->fetch(buffer.data());
            auto rext = tatami::consecutive_extractor<false>(&rhs, false, 0, rhs_col);
            size_t out_offset = static_cast<size_t>(r) * row_shift; // using offsets instead of directly adding to the pointer, to avoid forming an invalid address on the final iteration.

            for (RightIndex_ j = 0; j < rhs_col; ++j, out_offset += col_shift) {
                auto rptr = rext->fetch(rbuffer.data());
                output[out_offset] = std::inner_product(ptr, ptr + NC, rptr, static_cast<Output_>(0));
            }
        }

    }, NR, num_threads);
}

template<typename Value_, typename Index_, typename RightValue_, typename RightIndex_, typename Output_>
void dense_row_tatami_sparse(const tatami::Matrix<Value_, Index_>& matrix, const tatami::Matrix<RightValue_, RightIndex_>& rhs, Output_* output, size_t row_shift, size_t col_shift, int num_threads) {
    Index_ NR = matrix.nrow();
    Index_ NC = matrix.ncol();
    RightIndex_ rhs_col = rhs.ncol();

    tatami::parallelize([&](size_t, Index_ start, Index_ length) {
        auto ext = tatami::consecutive_extractor<false>(&matrix, true, start, length);
        std::vector<Value_> buffer(NC);
        std::vector<RightValue_> vbuffer(NC);
        std::vector<RightIndex_> ibuffer(NC);

        constexpr bool supports_specials = supports_special_values<Value_>();
        typename std::conditional<supports_specials, std::vector<Index_>, bool>::type specials;

        for (Index_ r = start, end = start + length; r < end; ++r) {
            auto ptr = ext->fetch(buffer.data());
            auto rext = tatami::consecutive_extractor<true>(&rhs, false, 0, rhs_col);
            size_t out_offset = static_cast<size_t>(r) * row_shift; // using offsets instead of directly adding to the pointer, to avoid forming an invalid address on the final iteration.

            if constexpr(supports_specials) {
                specials.clear();
                fill_special_index(NC, ptr, specials);
            }

            for (RightIndex_ j = 0; j < rhs_col; ++j, out_offset += col_shift) {
                auto range = rext->fetch(vbuffer.data(), ibuffer.data());

                if constexpr(supports_specials) {
                    if (specials.size()) {
                        output[out_offset] = special_dense_sparse_multiply<Output_>(specials, ptr, range);
                        continue;
                    }
                }

                output[out_offset] = dense_sparse_multiply<Output_>(ptr, range);
            }
        }

    }, NR, num_threads);
}

}

}

#endif
