#ifndef BYTEME_RAW_FILE_WRITER_HPP
#define BYTEME_RAW_FILE_WRITER_HPP

#include <vector>
#include <stdexcept>
#include <string>
#include <cstdio>
#include "Writer.hpp"
#include "SelfClosingFILE.hpp"

/**
 * @file RawFileWriter.hpp
 *
 * @brief Write bytes to a file without any extra transformations.
 */

namespace byteme {

/**
 * @brief Write bytes to a file.
 *
 * This class will write bytes to a file without any further transformations.
 * It is basically a simple wrapper around `FILE` structures, with correct closing and error checking.
 */
class RawFileWriter : public Writer {
public:
    /**
     * @param path Path to the file.
     * @param buffer_size Size of the buffer to use for writing.
     */
    RawFileWriter(const char* path, size_t buffer_size = 65536) : file(path, "wb") {
        if (std::setvbuf(file.handle, nullptr, _IOFBF, buffer_size)) {
            throw std::runtime_error("failed to set a buffer size for file writing");
        }
    }

    /**
     * @param path Path to the file.
     * @param buffer_size Size of the buffer to use for writing.
     */
    RawFileWriter(const std::string& path, size_t buffer_size = 65536) : RawFileWriter(path.c_str(), buffer_size) {}

public:
    using Writer::write;

    void write(const unsigned char* buffer, size_t n) {
        size_t ok = std::fwrite(buffer, sizeof(unsigned char), n, file.handle);
        if (ok < n) {
            throw std::runtime_error("failed to write raw binary file (fwrite error " + std::to_string(std::ferror(file.handle)) + ")");
        }
    }

    void finish() {
        if (std::fclose(file.handle)) {
            throw std::runtime_error("failed to close raw binary file");
        }
        file.handle = nullptr;
    }

private:
    SelfClosingFILE file;
};

}

#endif
