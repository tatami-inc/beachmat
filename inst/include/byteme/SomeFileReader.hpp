#ifndef BYTEME_SOME_FILE_READER_HPP
#define BYTEME_SOME_FILE_READER_HPP

#include "Reader.hpp"
#include "RawFileReader.hpp"
#include "GzipFileReader.hpp"
#include "magic_numbers.hpp"
#include <memory>
#include <cstdio>

/**
 * @file SomeFileReader.hpp
 *
 * @brief Read a possibly-Gzipped file.
 */

namespace byteme {

/**
 * @brief Read a file that may or may not be Gzipped.
 *
 * This class will automatically detect whether `path` refers to a text file or a Gzip-compressed file, based on its initial magic numbers.
 * After that, it will dispatch appropriately to `RawFileReader` or `GzipFileReader` respectively.
 */
class SomeFileReader : public Reader {
public:
    /**
     * @param path Path to the file.
     * @param buffer_size Size of the buffer to use for reading.
     */
    SomeFileReader(const char* path, size_t buffer_size = 65536) { 
        unsigned char header[3];
        size_t read;
        {
            SelfClosingFILE file(path, "rb");
            read = std::fread(header, sizeof(unsigned char), 3, file.handle);
        }

        if (is_gzip(header, read)) {
            source.reset(new GzipFileReader(path, buffer_size));
        } else {
            source.reset(new RawFileReader(path, buffer_size));
        }
    }

    /**
     * @param path Path to the file.
     * @param buffer_size Size of the buffer to use for reading.
     */
    SomeFileReader(const std::string& path, size_t buffer_size = 65536) : SomeFileReader(path.c_str(), buffer_size) {}

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
