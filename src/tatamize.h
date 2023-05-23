#ifndef TATAMIZE_H
#define TATAMIZE_H

// Order of includes is very important here.
#define TATAMI_R_PARALLELIZE_UNKNOWN
#include "tatami_r/parallelize.hpp"
#define TATAMI_CUSTOM_PARALLEL tatami_r::parallelize
#include "tatami_r/tatami_r.hpp"

#include "tatami/tatami.hpp"
#include <memory>
#include "Rcpp.h"

struct MatrixChan {
    MatrixChan(tatami::NumericMatrix* p) : ptr(p) {}
    MatrixChan(std::shared_ptr<tatami::NumericMatrix> p) : ptr(std::move(p)) {}
    std::shared_ptr<tatami::NumericMatrix> ptr;
    const tatami::NumericMatrix* get() const { return ptr.get(); } 
};

typedef Rcpp::XPtr<MatrixChan> MatrixChanPtr;

inline MatrixChanPtr new_MatrixChan(tatami::NumericMatrix* p) {
    return MatrixChanPtr(new MatrixChan(p), true); 
}

inline MatrixChanPtr new_MatrixChan(std::shared_ptr<tatami::NumericMatrix> p) {
    return MatrixChanPtr(new MatrixChan(std::move(p)), true); 
}

inline const tatami::NumericMatrix* extract_NumericMatrix(SEXP x) {
    MatrixChanPtr mat(x);
    return mat->get();
}

inline std::shared_ptr<tatami::NumericMatrix> extract_NumericMatrix_shared(SEXP x) {
    MatrixChanPtr mat(x);
    return mat->ptr;
}

#endif
