#ifndef BEACHMAT_DIM_CHECKER_H
#define BEACHMAT_DIM_CHECKER_H

#include "Rcpp.h"

/* Virtual base class for all reader/writer classes. */

namespace beachmat{

class dim_checker {
public:
    dim_checker() = default;
    dim_checker(size_t, size_t);

    virtual ~dim_checker() = default;
    dim_checker(const dim_checker&) = default;
    dim_checker& operator=(const dim_checker&) = default;
    dim_checker(dim_checker&&) = default;
    dim_checker& operator=(dim_checker&&) = default;

    size_t get_nrow() const;
    size_t get_ncol() const;

    // Helper functions that might be useful elsewhere.
    static void check_dimension(size_t, size_t, const char*);
    static void check_subset(size_t, size_t, size_t, const char*);
protected:
    size_t nrow=0, ncol=0;
    void fill_dims(const Rcpp::RObject&);

    void check_rowargs(size_t) const;
    void check_rowargs(size_t, size_t, size_t) const;

    void check_colargs(size_t) const;
    void check_colargs(size_t, size_t, size_t) const;

    void check_oneargs(size_t, size_t) const;

    void check_row_indices(Rcpp::IntegerVector::iterator, size_t);
    void check_col_indices(Rcpp::IntegerVector::iterator, size_t);
};

}

#endif
