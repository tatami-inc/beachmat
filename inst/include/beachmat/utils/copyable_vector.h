#ifndef BEACHMAT_COPYABLE_H
#define BEACHMAT_COPYABLE_H

#include "Rcpp.h"

namespace beachmat {

// Definition of a copyable vector, to use in classes.

template <class V> 
struct copyable_holder {
    V vec;

    copyable_holder(size_t n=0) : vec(n) {}
    ~copyable_holder() = default;
    
    // Ensuring that the copy generates a deep clone.
    copyable_holder(const copyable_holder& other) : vec(Rcpp::clone(other.vec)) {}
    copyable_holder& operator=(const copyable_holder& other) { 
        vec=Rcpp::clone(other.vec); 
        return *this;
    }

    // Using the copy constructor explicitly.
    copyable_holder(copyable_holder&& other) : vec(other.vec) {}
    copyable_holder& operator=(copyable_holder&& other) {
        vec=other.vec;
        return *this;
    }
};

}

#endif
