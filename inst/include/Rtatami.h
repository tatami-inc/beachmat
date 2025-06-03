#ifndef RTATAMI_H
#define RTATAMI_H

// Order of includes is very important here.
#define TATAMI_R_PARALLELIZE_UNKNOWN
#include "tatami_r/parallelize.hpp"
#define TATAMI_CUSTOM_PARALLEL tatami_r::parallelize
#include "tatami_r/tatami_r.hpp"

#include "tatami/tatami.hpp"
#include <memory>
#include "Rcpp.h"

namespace Rtatami {

/**
 * @brief Pointer to a **tatami** numeric matrix.
 *
 * The `tatami::NumericMatrix` is allowed to hold views on R-owned data, to avoid copies when moving from R to C++.
 * However, if garbage collection occurs on the R-owned data, the use of that data in C++ becomes invalid.
 * To avoid this, we hold the original R objects to protect them from R's garbage collector until this object is also destroyed.
 */
struct BoundNumericMatrix {
    /**
     * Pointer to a `tatami::NumericMatrix`.
     */
    std::shared_ptr<tatami::NumericMatrix> ptr;

    /**
     * The original R object.
     */
    Rcpp::RObject original;

    /**
     * @return Raw pointer to a `tatami::NumericMatrix`.
     */
    const tatami::NumericMatrix* get() const { return ptr.get(); } 
};

/**
 * A **Rcpp** external pointer to a `BoundNumericMatrix` object.
 */
typedef Rcpp::XPtr<BoundNumericMatrix> BoundNumericPointer;

/**
 * Create a new `BoundNumericMatrix` instance.
 * It is expected that functions will set `original` and `ptr` themselves before returning to the user.
 *
 * @return A `BoundNumericPointer` to a default-initialized (i.e., empty) `BoundNumericMatrix` object.
 */
inline BoundNumericPointer new_BoundNumericMatrix() {
    return Rcpp::XPtr<BoundNumericMatrix>(new BoundNumericMatrix, true); 
}

/**
 * Set or unset the parallel executor.
 * 
 * @param ptr An external pointer created by `beachmat::getExecutor()`, or NULL to unset the existing executor.
 */
inline void set_executor(SEXP ptr) {
    if (ptr == R_NilValue) {
        tatami_r::set_executor(NULL);
    } else {
        tatami_r::set_executor(Rcpp::XPtr<manticore::Executor>(ptr).get());
    }
} 

}

#endif
