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
protected:
    size_t nrow, ncol;
    void fill_dims(const Rcpp::RObject&);

    void check_rowargs(size_t) const;
    void check_rowargs(size_t, size_t, size_t) const;

    void check_colargs(size_t) const;
    void check_colargs(size_t, size_t, size_t) const;

    void check_oneargs(size_t, size_t) const;
};

}

#endif
