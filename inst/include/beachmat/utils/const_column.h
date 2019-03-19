#ifndef BEACHMAT_CONST_COLUMN_H
#define BEACHMAT_CONST_COLUMN_H

#include "../integer_matrix.h"
#include "../logical_matrix.h"
#include "../numeric_matrix.h"
#include "../character_matrix.h"

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
            raws=raw_structure<typename M::vector>(mat->get_nrow()); // repurposing the raw structure to hold some values.
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

    Rcpp::IntegerVector::iterator get_indices() const {
        if (Is_sparse) {
            return raws.get_structure_start();
        }
        if (indices.empty() || indices.back().size() < ptr->get_nrow()) {
            indices.push_back(Rcpp::IntegerVector(ptr->get_nrow()));
            std::iota(indices.back().begin(), indices.back().end(), 0);
        }
        return indices.back().begin()+last_start;
    }

private:
    M* ptr;
    raw_structure<typename M::vector> raws;
    bool Is_dense, Is_sparse;

    // Use a vector of vectors, to hold onto previous indexing vectors.
    // This avoids invalidating existing iterators if there are multiple 'mat's in play.
    static std::vector<Rcpp::IntegerVector> indices; 
    size_t last_start=0;
};

}

#endif
