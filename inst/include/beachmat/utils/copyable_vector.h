#ifndef BEACHMAT_COPYABLE_H
#define BEACHMAT_COPYABLE_H

#include "Rcpp.h"

namespace beachmat {

// Definition of a copyable vector, to use in classes.

template <class V> 
struct copyable_holder {
    copyable_holder(size_t n=0) : vec(n) {}
    ~copyable_holder() = default;
    copyable_holder(const copyable_holder& other) : vec(Rcpp::clone(other.vec)) {}
    copyable_holder& operator=(const copyable_holder& other) { 
        vec=Rcpp::clone(other.vec); 
        return *this;
    }
    copyable_holder(copyable_holder&&) = default;
    copyable_holder& operator=(copyable_holder&&) = default;
    V vec;
};

}

#endif
