#ifndef BEACHMAT_DIM_CHECKER_H
#define BEACHMAT_DIM_CHECKER_H

#include "Rcpp.h"
#include <stdexcept>

namespace beachmat {

class dim_checker {
public:
    dim_checker() {}

    virtual ~dim_checker() = default;
    dim_checker(const dim_checker&) = default;
    dim_checker& operator=(const dim_checker&) = default;
    dim_checker(dim_checker&&) = default;
    dim_checker& operator=(dim_checker&&) = default;

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

    size_t get_nrow() const { return nrow; }

    size_t get_ncol() const { return ncol; }
protected:
    size_t nrow=0, ncol=0;

    void fill_dims(Rcpp::RObject dims) {
        if (dims.sexp_type()!=INTSXP) {
            throw std::runtime_error("matrix dimensions should be an integer vector");
        }

        Rcpp::IntegerVector d(dims);
        if (d.size()!=2) {
            throw std::runtime_error("matrix dimensions should be of length 2");
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
