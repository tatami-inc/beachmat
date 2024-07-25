#ifndef TATAMI_MULT_SPARSE_COLUMN_HPP
#define TATAMI_MULT_SPARSE_COLUMN_HPP

#include <vector>
#include <cstdint>

#include "tatami/tatami.hpp"
#include "tatami_stats/tatami_stats.hpp"
#include "utils.hpp"

namespace tatami_mult {

namespace internal {

template<typename Value_, typename Index_, typename Output_>
void sparse_multiply_add(const tatami::SparseRange<Value_, Index_>& range, Index_ start, Output_ mult, Output_* optr) {
    for (Index_ r = 0; r < range.number; ++r) {
        optr[range.index[r] - start] += mult * range.value[r];
    }
}

template<typename Value_, typename Index_>
void expand_sparse_range(const tatami::SparseRange<Value_, Index_>& range, Index_ start, std::vector<Value_>& expanded) {
    for (Index_ k = 0; k < range.number; ++k) {
        expanded[range.index[k] - start] = range.value[k];
    }
}

template<typename Value_, typename Index_>
void reset_expanded_sparse_range(const tatami::SparseRange<Value_, Index_>& range, Index_ start, std::vector<Value_>& expanded) {
    for (Index_ k = 0; k < range.number; ++k) {
        expanded[range.index[k] - start] = 0;
    }
}

template<typename Value_, typename Index_, typename Right_, typename Output_>
void sparse_column_vector(const tatami::Matrix<Value_, Index_>& matrix, const Right_* rhs, Output_* output, int num_threads) {
    Index_ NR = matrix.nrow();
    Index_ NC = matrix.ncol();

    tatami::parallelize([&](size_t t, Index_ start, Index_ length) {
        auto ext = tatami::consecutive_extractor<true>(&matrix, false, 0, NC, start, length);
        std::vector<Value_> vbuffer(length);
        std::vector<Index_> ibuffer(length);

        tatami_stats::LocalOutputBuffer<Value_> store(t, start, length, output);
        auto optr = store.data();

        constexpr bool supports_specials = supports_special_values<Right_>();
        typename std::conditional<supports_specials, std::vector<Value_>, bool>::type expanded;

        for (Index_ c = 0; c < NC; ++c) {
            auto range = ext->fetch(vbuffer.data(), ibuffer.data());

            Output_ mult = rhs[c];
            if constexpr(supports_specials) {
                if (is_special(mult)) { 
                    // Expanding it to explicitly perform the multiplication with specials.
                    expanded.resize(length);
                    expand_sparse_range(range, start, expanded);

                    for (Index_ r = 0; r < length; ++r) {
                        optr[r] += expanded[r] * mult;
                    }

                    reset_expanded_sparse_range(range, start, expanded);
                    continue;
                }
            }

            sparse_multiply_add(range, start, mult, optr);
        }

        store.transfer();
    }, NR, num_threads);
}

template<typename Value_, typename Index_, typename Right_, typename Output_>
void sparse_column_vectors(const tatami::Matrix<Value_, Index_>& matrix, const std::vector<Right_*>& rhs, const std::vector<Output_*>& output, int num_threads) {
    Index_ NR = matrix.nrow();
    Index_ NC = matrix.ncol();
    size_t num_rhs = rhs.size();

    tatami::parallelize([&](size_t t, Index_ start, Index_ length) {
        auto ext = tatami::consecutive_extractor<true>(&matrix, false, 0, NC, start, length);
        std::vector<Value_> vbuffer(length);
        std::vector<Index_> ibuffer(length);
        auto stores = create_stores(t, start, length, output);

        constexpr bool supports_specials = supports_special_values<Right_>();
        typename std::conditional<supports_specials, std::vector<Value_>, bool>::type expanded;

        for (Index_ c = 0; c < NC; ++c) {
            auto range = ext->fetch(vbuffer.data(), ibuffer.data());
            bool has_expanded = false;

            for (size_t j = 0; j < num_rhs; ++j) {
                auto optr = stores[j].data();
                Output_ mult = rhs[j][c];

                if constexpr(supports_specials) {
                    if (is_special(mult)) { 
                        // Expanding it to explicitly perform the multiplication with specials.
                        if (!has_expanded) {
                            expanded.resize(length);
                            expand_sparse_range(range, start, expanded);
                            has_expanded = true;
                        }

                        for (Index_ r = 0; r < length; ++r) {
                            optr[r] += expanded[r] * mult;
                        }
                        continue;
                    }
                }

                sparse_multiply_add(range, start, mult, optr);
            }

            if constexpr(supports_specials) {
                if (has_expanded) {
                    reset_expanded_sparse_range(range, start, expanded);
                }
            }
        }
   
        for (auto& s : stores) {
            s.transfer();
        }
    }, NR, num_threads);
}

template<typename Value_, typename Index_, typename RightValue_, typename RightIndex_, typename Output_>
void sparse_column_tatami_dense(const tatami::Matrix<Value_, Index_>& matrix, const tatami::Matrix<RightValue_, RightIndex_>& rhs, Output_* output, size_t row_shift, size_t col_shift, int num_threads) {
    Index_ NR = matrix.nrow();
    Index_ NC = matrix.ncol();
    RightIndex_ rhs_col = rhs.ncol();

    tatami::parallelize([&](size_t t, Index_ start, Index_ length) {
        auto ext = tatami::consecutive_extractor<true>(&matrix, false, 0, NC, start, length);
        auto rext = tatami::consecutive_extractor<false>(&rhs, true, 0, NC); // remember, NC == rhs.nrow().
        std::vector<Value_> vbuffer(length);
        std::vector<Index_> ibuffer(length);
        std::vector<RightValue_> rbuffer(rhs_col);

        bool contiguous_output = (row_shift == 1);
        size_t mock_thread = (contiguous_output ? t : static_cast<size_t>(-1)); // avoid a direct right if it's not contiguous.
        auto stores = create_stores(mock_thread, start, length, output, rhs_col, col_shift);

        constexpr bool supports_specials = supports_special_values<RightValue_>();
        typename std::conditional<supports_specials, std::vector<Value_>, bool>::type expanded;

        for (Index_ c = 0; c < NC; ++c) {
            auto range = ext->fetch(vbuffer.data(), ibuffer.data());
            auto rptr = rext->fetch(rbuffer.data());
            bool has_expanded = false;

            for (RightIndex_ j = 0; j < rhs_col; ++j) {
                auto optr = stores[j].data();
                Output_ mult = rptr[j];

                if constexpr(supports_specials) {
                    if (is_special(mult)) { // expanding it to explicitly perform the multiplication with specials.
                        // Expanding it to explicitly perform the multiplication with specials.
                        if (!has_expanded) {
                            expanded.resize(length);
                            expand_sparse_range(range, start, expanded);
                            has_expanded = true;
                        }

                        for (Index_ r = 0; r < length; ++r) {
                            optr[r] += expanded[r] * mult;
                        }
                        continue;
                    }
                }

                sparse_multiply_add(range, start, mult, optr);
            }

            if constexpr(supports_specials) {
                if (has_expanded) {
                    reset_expanded_sparse_range(range, start, expanded);
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
void sparse_column_tatami_sparse(const tatami::Matrix<Value_, Index_>& matrix, const tatami::Matrix<RightValue_, RightIndex_>& rhs, Output_* output, size_t row_shift, size_t col_shift, int num_threads) {
    Index_ NR = matrix.nrow();
    Index_ NC = matrix.ncol();
    RightIndex_ rhs_col = rhs.ncol();

    tatami::parallelize([&](size_t t, Index_ start, Index_ length) {
        auto ext = tatami::consecutive_extractor<true>(&matrix, false, 0, NC, start, length);
        auto rext = tatami::consecutive_extractor<true>(&rhs, true, 0, NC); // remember, NC == rhs.nrow().
        std::vector<Value_> vbuffer(length);
        std::vector<Index_> ibuffer(length);
        std::vector<RightValue_> rvbuffer(rhs_col);
        std::vector<RightIndex_> ribuffer(rhs_col);

        bool contiguous_output = (row_shift == 1);
        size_t mock_thread = (contiguous_output ? t : static_cast<size_t>(-1)); // avoid a direct right if it's not contiguous.
        auto stores = create_stores(mock_thread, start, length, output, rhs_col, col_shift);

        // This time, we're checking for special values in the LHS because
        // we're potentially skipping its columns based on the sparsity of the
        // RHS and we need to know whether to NOT skip due to specials.
        constexpr bool supports_specials = supports_special_values<Value_>();
        typename std::conditional<supports_specials, std::vector<Index_>, bool>::type special_k;

        for (Index_ c = 0; c < NC; ++c) {
            auto range = ext->fetch(vbuffer.data(), ibuffer.data());
            auto rhs_range = rext->fetch(rvbuffer.data(), ribuffer.data());

            if constexpr(supports_specials) {
                special_k.clear();
                for (Index_ k = 0; k < range.number; ++k) {
                    if (is_special(range.value[k])) {
                        special_k.push_back(k);
                    }
                }
            }

            if constexpr(supports_specials) { // need separate multiplication to preserve the specials.
                if (special_k.size()) {
                    RightIndex_ rhs_k = 0; 
                    for (RightIndex_ j = 0; j < rhs_col; ++j) {
                        auto optr = stores[j].data();
                        if (rhs_k < rhs_range.number && j == rhs_range.index[rhs_k]) {
                            Output_ mult = rhs_range.value[rhs_k];
                            for (Index_ k = 0; k < range.number; ++k) {
                                optr[range.index[k] - start] += mult * range.value[k];
                            }
                            ++rhs_k;
                        } else {
                            for (auto k : special_k) {
                                optr[range.index[k] - start] += range.value[k] * static_cast<Output_>(0);
                            }
                        }
                    }
                    continue;
                }
            }

            for (RightIndex_ rhs_k = 0; rhs_k < rhs_range.number; ++rhs_k) {
                auto optr = stores[rhs_range.index[rhs_k]].data();
                Output_ mult = rhs_range.value[rhs_k];
                for (Index_ k = 0; k < range.number; ++k) {
                    optr[range.index[k] - start] += mult * range.value[k];
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
