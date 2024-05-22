#ifndef TATAMI_CHUNKED_TATAMI_CHUNKED_HPP
#define TATAMI_CHUNKED_TATAMI_CHUNKED_HPP

/**
 * @file tatami_chunked.hpp
 * @brief Umbrella header for chunked tatami matrices.
 */

#include "LruSlabCache.hpp"
#include "OracularSlabCache.hpp"
#include "OracularVariableSlabCache.hpp"
#include "OracularSubsettedSlabCache.hpp"

#include "SlabCacheStats.hpp"
#include "DenseSlabFactory.hpp"
#include "SparseSlabFactory.hpp"

#include "mock_dense_chunk.hpp"
#include "mock_sparse_chunk.hpp"
#include "CustomDenseChunkedMatrix.hpp"
#include "CustomSparseChunkedMatrix.hpp"

/**
 * @namespace tatami_chunked
 * @brief Methods to handle chunked tatami matrices.
 */
namespace tatami_chunked {}

#endif
