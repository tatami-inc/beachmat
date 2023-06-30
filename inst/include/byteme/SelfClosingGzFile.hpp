#ifndef BYTEME_SELF_CLOSING_GZFILE_HPP
#define BYTEME_SELF_CLOSING_GZFILE_HPP

#include <stdexcept>
#include <string>
#include "zlib.h"

namespace byteme {

struct SelfClosingGzFile {
    SelfClosingGzFile(const char* path, const char* mode) : handle(gzopen(path, mode)) {
        if (!handle) {
            throw std::runtime_error("failed to open file at '" + std::string(path) + "'");
        }
        return;
    }

    ~SelfClosingGzFile() {
        if (!closed) {
            gzclose(handle);
        }
        return;
    }

    SelfClosingGzFile(SelfClosingGzFile&& x) : handle(std::move(x.handle)) {
        x.closed = true;
    }

    SelfClosingGzFile& operator=(SelfClosingGzFile&& x) {
        handle = std::move(x.handle);
        x.closed = true;
        return *this;
    }

    // Delete the remaining constructors.
    SelfClosingGzFile(const SelfClosingGzFile&) = delete;
    SelfClosingGzFile& operator=(const SelfClosingGzFile&) = delete;

    bool closed = false;
    gzFile handle;
};

}

#endif
