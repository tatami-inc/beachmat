#ifndef BYTEME_CHUNKED_BUFFER_READER_HPP
#define BYTEME_CHUNKED_BUFFER_READER_HPP

#include <algorithm>
#include "Reader.hpp"

/**
 * @file ChunkedBufferReader.hpp
 *
 * @brief Read chunks of bytes from a raw buffer.
 */

namespace byteme {

/**
 * @brief Read chunks of bytes from a raw buffer.
 *
 * This is basically the same as `RawBufferReader` except that chunks of bytes are returned on every `load()` call.
 * It is primarily intended for use in tests of **byteme** callers, to ensure that downstream algorithms behave correctly with respect to chunked reads.
 */
struct ChunkedBufferReader : public byteme::Reader {
    /**
     * @param[in] buffer Pointer to an array of bytes, usually containing text.
     * @param length Length of the buffer.
     * @param chunk_size Size of each chunk in bytes.
     */
    ChunkedBufferReader(const unsigned char* buffer, size_t length, size_t chunk_size) : source(buffer), len(length), chunksize(chunk_size) {
        position = -chunksize;
    }

    /**
     * @param[in] buffer Pointer to an array of bytes, usually containing text.
     * @param length Length of the buffer.
     * @param chunk_size Size of each chunk in bytes.
     */
    ChunkedBufferReader(const char* buffer, size_t length, size_t chunk_size) : 
        ChunkedBufferReader(reinterpret_cast<const unsigned char*>(buffer), length, chunk_size) {}

public:
    bool load() {
        position += chunksize; 
        return (position < len);
    }

    const unsigned char * buffer () const {
        return source + position;
    }

    size_t available() const {
        return std::min(chunksize, len - position);
    }

private:
    const unsigned char* source;
    size_t len;
    size_t position;
    size_t chunksize;
};

}

#endif
