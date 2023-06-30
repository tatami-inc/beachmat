#ifndef BYTEME_GZIP_FILE_WRITER_HPP
#define BYTEME_GZIP_FILE_WRITER_HPP

#include <stdexcept>
#include <vector>
#include <string>
#include "Writer.hpp"

/**
 * @file GzipFileWriter.hpp
 *
 * @brief Write a Gzip-compressed file.
 */

namespace byteme {

/**
 * @brief Write uncompressed bytes to a Gzip-compressed file.
 *
 * This is basically a wrapper around Zlib's `gzFile` with correct closing and error checking.
 */
class GzipFileWriter : public Writer {
public:
    /**
     * @param path Path to the file.
     * @param compression_level Gzip compression level. 
     * @param buffer_size Size of the internal buffer. 
     */
    GzipFileWriter(const char* path, int compression_level = 6, size_t buffer_size = 65536) : gz(path, "wb") {
        if (gzbuffer(gz.handle, buffer_size)) {
            throw std::runtime_error("failed to set the Gzip compression buffer");
        }
        if (gzsetparams(gz.handle, compression_level, Z_DEFAULT_STRATEGY) != Z_OK) {
            throw std::runtime_error("failed to set the Gzip compression parameters");
        }
    }

    /**
     * @param path Path to the file.
     * @param compression_level Gzip compression level. 
     * @param buffer_size Size of the buffer to use for reading.
     */
    GzipFileWriter(const std::string& path, int compression_level = 6, size_t buffer_size = 65536) : GzipFileWriter(path.c_str(), compression_level, buffer_size) {}

public:
    using Writer::write;

    void write(const unsigned char* buffer, size_t n) {
        if (n) {
            size_t ok = gzwrite(gz.handle, buffer, n);
            if (ok != n) {
                throw std::runtime_error("failed to write to the Gzip-compressed file");
            }
        }
    }

    void finish() {
        gz.closed = true;
        if (gzclose(gz.handle) != Z_OK) {
            throw std::runtime_error("failed to close the Gzip-compressed file after writing");
        }
    }

private:
    SelfClosingGzFile gz;
};

}

#endif
