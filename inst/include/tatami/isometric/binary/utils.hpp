#ifndef TATAMI_ISOMETRIC_BINARY_HELPER_UTILS_H
#define TATAMI_ISOMETRIC_BINARY_HELPER_UTILS_H

namespace tatami {

template<bool must_have_both, typename Value_, typename Index_, class Function_>
Index_ delayed_binary_isometric_sparse_operation(
    const SparseRange<Value_, Index_>& left,
    const SparseRange<Value_, Index_>& right, 
    Value_* value_buffer, 
    Index_* index_buffer, 
    bool needs_value, 
    bool needs_index, 
    Function_ fun)
{
    Index_ lcount = 0, rcount = 0, output = 0;

    auto advance_left = [&]() -> void {
        if (needs_value) {
            value_buffer[output] = left.value[lcount];
            fun(value_buffer[output], 0);
        }
        if (needs_index) {
            index_buffer[output] = left.index[lcount];
        }
        ++output;
        ++lcount;
    };

    auto advance_right = [&]() -> void {
        if (needs_value) {
            value_buffer[output] = 0;
            fun(value_buffer[output], right.value[rcount]);
        }
        if (needs_index) {
            index_buffer[output] = right.index[rcount];
        }
        ++rcount;
        ++output;
    };

    while (lcount < left.number && rcount < right.number) {
        if (left.index[lcount] < right.index[rcount]) {
            if constexpr(!must_have_both) {
                advance_left();
            } else {
                ++lcount;
            }

        } else if (left.index[lcount] > right.index[rcount]) {
            if constexpr(!must_have_both) {
                advance_right();
            } else {
                ++rcount;
            }

        } else {
            if (needs_value) {
                value_buffer[output] = left.value[lcount];
                fun(value_buffer[output], right.value[rcount]);
            }
            if (needs_index) {
                index_buffer[output] = right.index[rcount];
            }
            ++lcount;
            ++rcount;
            ++output;
        }
    }

    if constexpr(!must_have_both) {
        while (lcount < left.number) {
            advance_left();
        }

        while (rcount < right.number) {
            advance_right();
        }
    }

    return output;
}

}

#endif
