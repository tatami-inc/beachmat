#ifndef BYTEME_BYTEME_HPP
#define BYTEME_BYTEME_HPP

#include "Reader.hpp"
#include "RawBufferReader.hpp"
#include "RawFileReader.hpp"
#include "IstreamReader.hpp"
#include "ChunkedBufferReader.hpp"

#include "Writer.hpp"
#include "RawBufferWriter.hpp"
#include "RawFileWriter.hpp"
#include "OstreamWriter.hpp"

#include "PerByte.hpp"

#if __has_include("zlib.h")
#include "GzipFileReader.hpp"
#include "ZlibBufferReader.hpp"

#include "SomeBufferReader.hpp"
#include "SomeFileReader.hpp"

#include "GzipFileWriter.hpp"
#include "ZlibBufferWriter.hpp"
#endif

/**
 * @file byteme.hpp
 * @brief Umbrella header for all **byteme** classes
 *
 * If ZLib is not available, all of the Zlib-related headers are omitted.
 * This will skip classes such as the `GzipFileReader` and `SomeBufferReader`.
 */

/**
 * @namespace byteme
 * @brief Simple byte readers and writers
 */
namespace byteme {}

#endif
