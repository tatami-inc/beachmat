#ifndef BEACHMAT_CONST_COLUMN_H
#define BEACHMAT_CONST_COLUMN_H

#include "Rcpp.h"
#include "raw_structure.h"
#include <algorithm>

namespace beachmat {

/* A convenience class to obtain constant columns from a given matrix,
 * while supporting native sparse representations more-or-less transparently.
 */

template<class M>
class const_column {
public:
    const_column(M* mat, bool allow_sparse=true) : ptr(mat), raws(mat->set_up_raw()),
        Is_dense(mat->col_raw_type()=="dense"), Is_sparse(allow_sparse && mat->col_raw_type()=="sparse") 
    {
        if (!Is_dense && !Is_sparse) {
            // repurposing the raw structure to hold some values.
            raws=raw_structure<typename M::vector>(true, false, mat->get_nrow(), 0); 
        }
        return;
    }

    bool is_sparse () const { return Is_sparse; }

    bool is_dense () const { return Is_dense; }

    void fill(size_t c, size_t first, size_t last) {
        if (Is_dense || Is_sparse) {
            ptr->get_col_raw(c, raws, first, last);
        } else {
            ptr->get_col(c, raws.get_values_start(), first, last);
        }
        if (!Is_sparse) {
            last_start=first;
        }
        return;
    }

    void fill(size_t c) {
        fill(c, 0, ptr->get_nrow());
        return;
    }

    size_t get_n () const {
        return (Is_sparse ? raws.get_n() : ptr->get_nrow());
    }

    typename M::vector::iterator get_values() const {
        return raws.get_values_start();
    }

    void fill_indices() {
        if (!indices.size()) {
            indices=Rcpp::IntegerVector(ptr->get_nrow());
            std::iota(indices.begin(), indices.end(), 0);
        }
        return;
    }

    Rcpp::IntegerVector::iterator get_indices() {
        if (Is_sparse) {
            return raws.get_structure_start();
        }
        fill_indices();
        return indices.begin()+last_start;
    }
private:
    M* ptr;
    raw_structure<typename M::vector> raws;
    bool Is_dense, Is_sparse;

    Rcpp::IntegerVector indices; // deliberately copyable; values won't change, and reassignment won't matter.
    size_t last_start=0;
};

}

#endif
