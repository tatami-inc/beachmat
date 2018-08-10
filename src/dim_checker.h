#ifndef BEACHMAT_DIM_CHECKER_H
#define BEACHMAT_DIM_CHECKER_H

#include "beachmat.h"

/* Virtual base class for all reader/writer classes. */

namespace beachmat{

class dim_checker {
public:
    dim_checker();
    dim_checker(size_t, size_t);
    virtual ~dim_checker();
    size_t get_nrow() const;
    size_t get_ncol() const;

    // Helper functions that might be useful elsewhere.
    static void check_dimension(size_t, size_t, const char*);
    static void check_subset(size_t, size_t, size_t, const char*);
protected:
    size_t nrow, ncol;
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
