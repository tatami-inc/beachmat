#ifndef BYTEME_RAW_BUFFER_READER_HPP
#define BYTEME_RAW_BUFFER_READER_HPP

#include <algorithm>
#include "Reader.hpp"

/**
 * @file RawBufferReader.hpp
 *
 * @brief Read bytes from a raw buffer without any extra transformations.
 */

namespace byteme {

/**
 * @brief Read bytes from a raw buffer, usually text.
 *
 * This is a wrapper around an input buffer, provided for consistency with the other `*Reader` classes.
 * We assume that the lifetime of the data in the `buffer` pointer exceeds the lifetime of the instance.
 */
class RawBufferReader : public Reader {
public:
    /**
     * @param[in] buffer Pointer to an array of bytes, usually containing text.
     * @param length Length of the buffer.
     */
    RawBufferReader(const unsigned char* buffer, size_t length) : buffer_(buffer), len_(length) {}

    /**
     * @param[in] buffer Pointer to an array of bytes, usually containing text.
     * @param length Length of the buffer.
     */
    RawBufferReader(const char* buffer, size_t length) : buffer_(reinterpret_cast<const unsigned char*>(buffer)), len_(length) {}

public:
    bool load() {
        if (used) {
            return false;
        }
        used = true;
        return true;
    }

    const unsigned char* buffer() const {
        return buffer_;
    }

    size_t available() const {
        return len_;
    }

private:
    const unsigned char* buffer_;
    size_t len_;
    bool used = false;
};

}

#endif
