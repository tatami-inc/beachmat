#ifndef BYTEME_RAW_BUFFER_WRITER_HPP
#define BYTEME_RAW_BUFFER_WRITER_HPP

#include "Writer.hpp"
#include <vector>

/**
 * @file RawBufferWriter.hpp
 *
 * @brief Write bytes to a raw buffer without any extra transformations.
 */

namespace byteme {

/**
 * @brief Write bytes to a raw buffer.
 *
 * This class will append bytes to an internal instance of a `std::vector` without any further transformations.
 * Not much else to say here.
 */
class RawBufferWriter : public Writer {
public:
    /**
     * @param n Initial size of the output buffer to reserve.
     */
    RawBufferWriter(size_t n = 0) {
        output.reserve(n);
    }

    using Writer::write;

    void write(const unsigned char* buffer, size_t n) {
        output.insert(output.end(), buffer, buffer + n);
    }

    void finish() {}

    /**
     * Contents of the output buffer.
     * This should only be accessed after `finish()` is called.
     */
    std::vector<unsigned char> output;
};

}

#endif
