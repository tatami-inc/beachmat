#ifndef TATAMI_MULT_DENSE_COLUMN_HPP
#define TATAMI_MULT_DENSE_COLUMN_HPP

#include <vector>

#include "tatami/tatami.hpp"
#include "tatami_stats/tatami_stats.hpp"
#include "utils.hpp"

namespace tatami_mult {

namespace internal {

template<typename Value_, typename Index_, typename Right_, typename Output_>
void dense_column_vector(const tatami::Matrix<Value_, Index_>& matrix, const Right_* rhs, Output_* output, int num_threads) {
    Index_ NR = matrix.nrow();
    Index_ NC = matrix.ncol();

    tatami::parallelize([&](size_t t, Index_ start, Index_ length) {
        auto ext = tatami::consecutive_extractor<false>(&matrix, false, 0, NC, start, length);
        std::vector<Value_> buffer(NC);
        tatami_stats::LocalOutputBuffer<Value_> store(t, start, length, output);
        auto optr = store.data();

        for (Index_ c = 0; c < NC; ++c) {
            auto ptr = ext->fetch(buffer.data());
            Output_ mult = rhs[c];
            for (Index_ r = 0; r < length; ++r) {
                optr[r] += mult * ptr[r];
            }
        }

        store.transfer();
    }, NR, num_threads);
}

template<typename Value_, typename Index_, typename Right_, typename Output_>
void dense_column_vectors(const tatami::Matrix<Value_, Index_>& matrix, const std::vector<Right_*>& rhs, const std::vector<Output_*>& output, int num_threads) {
    Index_ NR = matrix.nrow();
    Index_ NC = matrix.ncol();
    size_t num_rhs = rhs.size();

    tatami::parallelize([&](size_t t, Index_ start, Index_ length) {
        auto ext = tatami::consecutive_extractor<false>(&matrix, false, 0, NC, start, length);
        std::vector<Value_> buffer(length);
        auto stores = create_stores(t, start, length, output);

        for (Index_ c = 0; c < NC; ++c) {
            auto ptr = ext->fetch(buffer.data());

            for (size_t j = 0; j < num_rhs; ++j) {
                auto optr = stores[j].data();
                Output_ mult = rhs[j][c];
                for (Index_ r = 0; r < length; ++r) {
                    optr[r] += mult * ptr[r];
                }
            }
        }
   
        for (auto& s : stores) {
            s.transfer();
        }
    }, NR, num_threads);
}

template<typename Value_, typename Index_, typename RightValue_, typename RightIndex_, typename Output_>
void dense_column_tatami_dense(const tatami::Matrix<Value_, Index_>& matrix, const tatami::Matrix<RightValue_, RightIndex_>& rhs, Output_* output, size_t row_shift, size_t col_shift, int num_threads) {
    Index_ NR = matrix.nrow();
    Index_ NC = matrix.ncol();
    RightIndex_ rhs_col = rhs.ncol();

    tatami::parallelize([&](size_t t, Index_ start, Index_ length) {
        auto ext = tatami::consecutive_extractor<false>(&matrix, false, 0, NC, start, length);
        auto rext = tatami::consecutive_extractor<false>(&rhs, true, 0, NC); // remember, NC == rhs.nrow().
        std::vector<Value_> buffer(length);
        std::vector<RightValue_> rbuffer(rhs_col);

        bool contiguous_output = (row_shift == 1);
        size_t mock_thread = (contiguous_output ? t : static_cast<size_t>(-1)); // avoid a direct right if it's not contiguous.
        auto stores = create_stores(mock_thread, start, length, output, rhs_col, col_shift);

        for (Index_ c = 0; c < NC; ++c) {
            auto ptr = ext->fetch(buffer.data());
            auto rptr = rext->fetch(rbuffer.data());

            for (RightIndex_ j = 0; j < rhs_col; ++j) {
                auto optr = stores[j].data();
                Output_ mult = rptr[j];
                for (Index_ r = 0; r < length; ++r) {
                    optr[r] += mult * ptr[r];
                }
            }
        }

        if (contiguous_output) {
            for (auto& s : stores) {
                s.transfer();
            }
        } else {
            non_contiguous_transfer(stores, start, length, output, row_shift, col_shift);
        }
    }, NR, num_threads);
}

template<typename Value_, typename Index_, typename RightValue_, typename RightIndex_, typename Output_>
void dense_column_tatami_sparse(const tatami::Matrix<Value_, Index_>& matrix, const tatami::Matrix<RightValue_, RightIndex_>& rhs, Output_* output, size_t row_shift, size_t col_shift, int num_threads) {
    Index_ NR = matrix.nrow();
    Index_ NC = matrix.ncol();
    RightIndex_ rhs_col = rhs.ncol();

    tatami::parallelize([&](size_t t, Index_ start, Index_ length) {
        auto ext = tatami::consecutive_extractor<false>(&matrix, false, 0, NC, start, length);
        auto rext = tatami::consecutive_extractor<true>(&rhs, true, 0, NC); // remember, NC == rhs.nrow().
        std::vector<Value_> buffer(length);
        std::vector<RightValue_> vbuffer(rhs_col);
        std::vector<RightIndex_> ibuffer(rhs_col);

        bool contiguous_output = (row_shift == 1);
        size_t mock_thread = (contiguous_output ? t : static_cast<size_t>(-1)); // avoid a direct right if it's not contiguous.
        auto stores = create_stores(mock_thread, start, length, output, rhs_col, col_shift);

        constexpr bool supports_specials = supports_special_values<Value_>();
        typename std::conditional<supports_specials, std::vector<Index_>, bool>::type specials;

        for (Index_ c = 0; c < NC; ++c) {
            auto ptr = ext->fetch(buffer.data());
            auto range = rext->fetch(vbuffer.data(), ibuffer.data());

            if constexpr(supports_specials) { // need separate multiplication to preserve the specials.
                specials.clear();
                fill_special_index(length, ptr, specials);

                if (specials.size()) {
                    RightIndex_ k = 0; 
                    for (RightIndex_ j = 0; j < rhs_col; ++j) {
                        auto optr = stores[j].data();
                        if (k < range.number && j == range.index[k]) {
                            Output_ mult = range.value[k];
                            for (Index_ r = 0; r < length; ++r) {
                                optr[r] += mult * ptr[r];
                            }
                            ++k;
                        } else {
                            for (auto s : specials) {
                                optr[s] += ptr[s] * static_cast<Output_>(0);
                            }
                        }
                    }
                    continue;
                }
            }

            for (RightIndex_ k = 0; k < range.number; ++k) {
                auto optr = stores[range.index[k]].data();
                Output_ mult = range.value[k];
                for (Index_ r = 0; r < length; ++r) {
                    optr[r] += mult * ptr[r];
                }
            }
        }

        if (contiguous_output) {
            for (auto& s : stores) {
                s.transfer();
            }
        } else {
            non_contiguous_transfer(stores, start, length, output, row_shift, col_shift);
        }
    }, NR, num_threads);
}

}

}

#endif
