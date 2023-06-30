#ifndef BYTEME_SOME_BUFFER_READER_HPP
#define BYTEME_SOME_BUFFER_READER_HPP

#include "Reader.hpp"
#include "RawBufferReader.hpp"
#include "ZlibBufferReader.hpp"
#include "magic_numbers.hpp"
#include <memory>
#include <cstdio>

/**
 * @file SomeBufferReader.hpp
 *
 * @brief Read a possibly-Gzipped or Zlibbed buffer.
 */

namespace byteme {

/**
 * @brief Read a buffer that may or may not be Gzip/Zlib-compressed.
 *
 * This class will automatically detect whether `buffer` refers to a text or Gzip/Zlib-compressed buffer, based on the initial magic numbers.
 * After that, it will dispatch appropriately to `RawBufferReader` or `ZlibBufferReader` respectively.
 */
class SomeBufferReader : public Reader {
public:
    /**
     * @param[in] buffer Pointer to an array containing the possibly compressed data.
     * @param len Length of the `buffer` array.
     * @param buffer_size Size of the buffer to use for decompression.
     */
    SomeBufferReader(const unsigned char* buffer, size_t len, size_t buffer_size = 65536) {
        if (is_zlib(buffer, len) || is_gzip(buffer, len)) {
            source.reset(new ZlibBufferReader(buffer, len, 3, buffer_size));
        } else {
            source.reset(new RawBufferReader(buffer, len));
        }
    }

    /**
     * @param[in] buffer Pointer to an array containing the possibly compressed data.
     * @param len Length of the `buffer` array.
     * @param buffer_size Size of the buffer to use for decompression.
     */
    SomeBufferReader(const char* buffer, size_t len, size_t buffer_size = 65536) :
        SomeBufferReader(reinterpret_cast<const unsigned char*>(buffer), len, buffer_size) {}

public:
    bool load() {
        return source->load();
    }

    const unsigned char* buffer() const {
        return source->buffer();
    }

    size_t available() const {
        return source->available();
    }

private:
    std::unique_ptr<Reader> source;
};

}

#endif
