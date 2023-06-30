#ifndef BYTEME_RAW_FILE_READER_HPP
#define BYTEME_RAW_FILE_READER_HPP

#include <vector>
#include <stdexcept>
#include <string>
#include <cstdio>
#include "Reader.hpp"
#include "SelfClosingFILE.hpp"

/**
 * @file RawFileReader.hpp
 *
 * @brief Read a file without any extra transformations.
 */

namespace byteme {

/**
 * @brief Read bytes from a file, usually text.
 *
 * This is basically a simple wrapper around `FILE` structures, with correct closing and error checking.
 * Mostly provided because I always forget how to interact with `ifstream` objects when I want a sequence of bytes.
 */
class RawFileReader : public Reader {
private:
    /**
     * @cond
     */
    friend class SomeFileReader;
    /**
     * @endcond
     */

public:
    /**
     * @param path Path to the file.
     * @param buffer_size Size of the buffer to use for reading.
     */
    RawFileReader(const char* path, size_t buffer_size = 65536) : file(path, "rb"), buffer_(buffer_size) {}

    /**
     * @param path Path to the file.
     * @param buffer_size Size of the buffer to use for reading.
     */
    RawFileReader(const std::string& path, size_t buffer_size = 65536) : RawFileReader(path.c_str(), buffer_size) {}

    bool load() {
        if (!okay) {
            return false;
        }

        auto& handle = file.handle;
        read = std::fread(buffer_.data(), sizeof(unsigned char), buffer_.size(), handle);

        if (read < buffer_.size()) {
            if (std::feof(handle)) {
                okay = false;
            } else {
                throw std::runtime_error("failed to read raw binary file (fread error " + std::to_string(std::ferror(handle)) + ")");
            }
        }

        return true;
    }

    const unsigned char* buffer() const {
        return buffer_.data();
    }

    size_t available() const {
        return read;
    }

private:
    SelfClosingFILE file;
    std::vector<unsigned char> buffer_;
    size_t read = 0;
    bool okay = true;
};

}

#endif
