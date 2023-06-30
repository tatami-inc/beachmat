#ifndef BYTEME_READER_HPP
#define BYTEME_READER_HPP

/**
 * @file Reader.hpp
 *
 * @brief Read an input source.
 */

namespace byteme {

/**
 * @brief Virtual class for reading bytes from a source.
 */
class Reader {
public:
    virtual ~Reader() = default;

    /**
     * Read the next chunk of bytes from the input source. 
     * To read the entire source, this function should be called repeatedly until `false` is returned.
     *
     * @return Boolean indicating whether the read was successful.
     * If `false`, it can be assumed that the end of the source was reached.
     */
    virtual bool load() = 0;

    /**
     * This method should only be called after `load()` has been called and returns `true`.
     *
     * @return Pointer to the start of an array containing the available bytes.
     * The number of available bytes is provided in `available()`.
     */
    virtual const unsigned char* buffer() const = 0;

    /**
     * This method should only be called after `load()` has been called and returns `true`.
     * The return value is generally expected to be positive; however, it is possible to return a zero.
     * Note that zero values should not be interpreted as the end of the source, which is strictly only defined by `load()` returning `false`.
     *
     * @return Number of available bytes in `buffer()`.
     */
    virtual size_t available() const = 0;
};

}

#endif
