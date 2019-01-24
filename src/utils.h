#ifndef BEACHMAT_UTILS_H
#define BEACHMAT_UTILS_H

#include "Rcpp.h"

#include <string>
#include <sstream>
#include <utility>
#include <tuple>

namespace beachmat { 

/* String-related helper functions */

std::string make_to_string(const Rcpp::RObject&);

template <class L, class R>
std::string combine_strings(const L& left, const R& right) {
    std::stringstream err;
    err << left << right;
    return err.str();    
}

void throw_custom_error(const std::string&, const std::string&, const std::string&);

// Class checking.

std::string get_class(const Rcpp::RObject&);

std::pair<std::string, std::string> get_class_package(const Rcpp::RObject&);

Rcpp::RObject get_safe_slot(const Rcpp::RObject&, const std::string&);

std::string check_Matrix_class (const Rcpp::RObject&, const std::string&);

// Type checking.

std::string translate_type(int);

int find_sexp_type (const Rcpp::RObject&);

// Dispatch function for external access check.

bool has_external_support (const Rcpp::RObject&);

// Matrix type enumeration.

enum matrix_type { SIMPLE, HDF5, SPARSE, DENSE, DELAYED, UNKNOWN, EXTERNAL };

// Typedef for the indexing tuple.

template<class V>
using const_col_indexed_info=std::tuple<size_t, Rcpp::IntegerVector::iterator, typename V::iterator>;

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
