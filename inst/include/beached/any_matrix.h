#ifndef BEACHMAT_ANY_MATRIX_H
#define BEACHMAT_ANY_MATRIX_H

#include "Rcpp.h"
#include <stdexcept>

namespace beachmat {

template <class T>
class any_matrix {
public:
    any_matrix() {}

    virtual ~any_matrix() = default;
    any_matrix(const any_matrix&) = default;
    any_matrix& operator=(const any_matrix&) = default;
    any_matrix(any_matrix&&) = default;
    any_matrix& operator=(any_matrix&&) = default;

    virtual T* get_col(size_t, T*, size_t, size_t);
    virtual T* get_row(size_t, T*, size_t, size_t);

    T* get_col(size_t r, T* work) {
        return get_col(r, work, 0, this->nrow);
    }

    T* get_row(size_t c, T* work) {
        return get_row(c, work, 0, this->ncol);        
    }

    // Helper functions that might be useful elsewhere.
    static void check_dimension(size_t i, size_t dim, const std::string& msg) {
        if (i >= dim) {
            throw std::runtime_error(msg + " index out of range");
        }
        return;
    }

    static void check_subset(size_t first, size_t last, size_t dim, const std::string& msg) {
         if (last < first) {
            throw std::runtime_error(msg + " start index is greater than " + msg + " end index");
         } else if (last > dim) {
             throw std::runtime_error(msg + " end index out of range");
         }
         
         return;    
    }
protected:
    size_t nrow=0, ncol=0;

    void fill_dims(const Rcpp::RObject& dims) {
        Rcpp::IntegerVector d;
        if (dims.sexp_type()!=d.sexp_type() || (d=dims).size()!=2) {
            throw std::runtime_error("matrix dimensions should be an integer vector of length 2");
        }
        if (d[0]<0 || d[1]<0) {
            throw std::runtime_error("dimensions should be non-negative");
        }
        nrow=d[0];
        ncol=d[1];
        return;
    }

    void check_rowargs(size_t r) const {
        dim_checker::check_dimension(r, nrow, "row");
        return;
    }

    void check_rowargs(size_t r, size_t first, size_t last) const {
        check_rowargs(r);
        dim_checker::check_subset(first, last, ncol, "column");
        return;
    }

    void check_colargs(size_t c) const {
        dim_checker::check_dimension(c, ncol, "column");
        return;
    }

    void check_colargs(size_t c, size_t first, size_t last) const {
        check_colargs(c);
        dim_checker::check_subset(first, last, nrow, "row");
        return;
    }
};

}

#endif
