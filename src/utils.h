#ifndef BEACHMAT_UTILS_H
#define BEACHMAT_UTILS_H

#include "beachmat.h"

namespace beachmat { 

std::string make_to_string(const Rcpp::RObject&);

std::string make_to_string(const char*, const char*);

void throw_custom_error(const std::string&, const std::string&, const std::string&);

// Class checking.

std::string get_class(const Rcpp::RObject&);

Rcpp::RObject get_safe_slot(const Rcpp::RObject&, const std::string&);

std::string check_Matrix_class (const Rcpp::RObject&, const std::string&);

// Type checking.

std::string translate_type(int);

int find_sexp_type (const Rcpp::RObject&);

// Matrix type enumeration.

enum matrix_type { SIMPLE, HDF5, SPARSE, DENSE, DELAYED, UNKNOWN };

// Typedef for the indexing tuple.

template<class V>
using const_col_indexed_info=std::tuple<size_t, Rcpp::IntegerVector::iterator, typename V::iterator>;

// Definition of a copyable vector, to use in classes.

template <class V> 
struct copyable_holder {
    copyable_holder(size_t n=0) : vec(n) {}
    ~copyable_holder() {};
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
