#include "delayed_matrix.h"
#include "unknown_matrix.h"

namespace beachmat {

/* Implementing methods for the 'delayed_coord_transformer' class */

template<typename T, class V>
delayed_coord_transformer<T, V>::delayed_coord_transformer() : transposed(false), byrow(false), bycol(false), 
        delayed_nrow(0), delayed_ncol(0) {}

template<typename T, class V>
template<class M>
delayed_coord_transformer<T, V>::delayed_coord_transformer(M mat) : transposed(false), byrow(false), bycol(false), 
        delayed_nrow(mat->get_nrow()), delayed_ncol(mat->get_ncol()) {}

template<typename T, class V>
template<class M>
delayed_coord_transformer<T, V>::delayed_coord_transformer(const Rcpp::RObject& in, M mat) : transposed(false), byrow(false), bycol(false),
        delayed_nrow(mat->get_nrow()), delayed_ncol(mat->get_ncol()), tmp(std::max(delayed_nrow, delayed_ncol)) {
   
    check_DelayedMatrix(in);
    if (!only_delayed_coord_changes(in)) { 
        throw std::runtime_error("'delayed_coord_transformer' called with non-empty delayed operations");
    }
    
    Rcpp::List indices(get_safe_slot(in, "index"));
    const size_t original_nrow(mat->get_nrow()), original_ncol(mat->get_ncol()); 

    // Checking indices for rows.
    Rcpp::RObject rowdex(indices[0]);
    byrow=!rowdex.isNULL();
    if (byrow){ 
        if (rowdex.sexp_type()!=INTSXP) {
            throw std::runtime_error("index vector should be integer");
        }

        Rcpp::IntegerVector rx(rowdex);
        row_index.insert(row_index.end(), rx.begin(), rx.end());
        for (auto& r : row_index) { 
            --r; // 0-based indices.
        }
        delayed_nrow=row_index.size();

        // If the indices are all consecutive from 0 to N, we turn byrow=false.
        if (delayed_nrow && row_index.front()==0 && delayed_nrow==original_nrow) {
            int count=0;
            byrow=false;
            for (auto r : row_index) {
                if (r!=count) {
                    byrow=true;
                    break;
                }
                ++count;
            }
        } 
    }

    // Checking indices for columns.
    Rcpp::RObject coldex(indices[1]);
    bycol=!coldex.isNULL();
    if (bycol){ 
        if (coldex.sexp_type()!=INTSXP) {
            throw std::runtime_error("index vector should be integer");
        }

        Rcpp::IntegerVector cx(coldex);
        col_index.insert(col_index.end(), cx.begin(), cx.end());
        for (auto& c : col_index) { 
            --c; // 0-based indices.
        }
        delayed_ncol=col_index.size();

        // If the indices are all consecutive from 0 to N, we turn bycol=false.
        if (delayed_ncol && col_index.front()==0 && delayed_ncol==original_ncol) {
            int count=0;
            bycol=false;
            for (auto c : col_index) {
                if (c!=count) {
                    bycol=true;
                    break;
                }
                ++count;
            }
        } 
    }
    
    // Checking transposition by peering into the "SeedDimChecker" class.
    Rcpp::RObject seed=get_safe_slot(in, "seed");
    if (seed.isS4() && get_class(seed)=="SeedDimPicker") { 
        Rcpp::IntegerVector dimorder(get_safe_slot(seed, "dim_combination"));
        if (dimorder.size()!=2) {
            throw std::runtime_error("'dim_combination' should be an integer vector of length 2");
        }
        transposed=(dimorder[0]==2);

        // As the row/column indices refer to the matrix BEFORE transposition.
        if (transposed) {
            std::swap(delayed_nrow, delayed_ncol);
        }
    }

    return;
}

template<typename T, class V>
size_t delayed_coord_transformer<T, V>::get_nrow() const{ 
    return delayed_nrow;
}

template<typename T, class V>
size_t delayed_coord_transformer<T, V>::get_ncol() const{ 
    return delayed_ncol;
}

template<typename T, class V>
template<class M, class Iter>
void delayed_coord_transformer<T, V>::get_row(M mat, size_t r, Iter out, size_t first, size_t last) {
    if (transposed) {
        r=transform_col(r);

        // Column extraction, first/last refer to rows.
        if (byrow) {
            mat->get_col(r, tmp.vec.begin());
            reallocate_col(first, last, out);
        } else {
            mat->get_col(r, out, first, last);
        }
    } else {
        r=transform_row(r);

        // Row extraction, first/last refer to columns.
        if (bycol) {
            mat->get_row(r, tmp.vec.begin());
            reallocate_row(first, last, out);
        } else {
            mat->get_row(r, out, first, last);
        }
    }
    return;
}

template<typename T, class V>
template<class M, class Iter>
void delayed_coord_transformer<T, V>::get_col(M mat, size_t c, Iter out, size_t first, size_t last) {
    if (transposed) {
        c=transform_row(c);

        // Row extraction, first/last refer to columns.
        if (bycol) {
            mat->get_row(c, tmp.vec.begin());
            reallocate_row(first, last, out);
        } else {
            mat->get_row(c, out, first, last);
        }
    } else {
        c=transform_col(c);

        // Column extraction, first/last refer to rows.
        if (byrow) {
            mat->get_col(c, tmp.vec.begin());
            reallocate_col(first, last, out);
        } else {
            mat->get_col(c, out, first, last);
        }
    }
    return;
}

template<typename T, class V>
template<class M>
T delayed_coord_transformer<T, V>::get(M mat, size_t r, size_t c) {
    if (transposed) {
        return mat->get(transform_row(c), transform_col(r));
    } else {
        return mat->get(transform_row(r), transform_col(c));
    }
}


template<typename T, class V>
size_t delayed_coord_transformer<T, V>::transform_row(size_t r) const {
    if (byrow) {
        if (r<0 || r>=row_index.size()) { 
            throw std::runtime_error("row indices out of range");
        }
        r=row_index[r];
    }
    return r;
}

template<typename T, class V>
size_t delayed_coord_transformer<T, V>::transform_col(size_t c) const {
    if (bycol) {
        if (c<0 || c>=col_index.size()) { 
            throw std::runtime_error("column indices out of range");
        }
        c=col_index[c];
    }
    return c;
}

template<typename T, class V>
template<class Iter>
void delayed_coord_transformer<T, V>::reallocate_row(size_t first, size_t last, Iter out) {
    if (first>last || last>col_index.size()) { 
        throw std::runtime_error("invalid column start/end indices");
    }
    
    const V& holding=tmp.vec;
    auto cIt=col_index.begin()+first, end=col_index.begin()+last;
    while (cIt!=end) {
        (*out)=holding[*cIt];
        ++out;
        ++cIt;
    }
    return;
}

template<typename T, class V>
template<class Iter>
void delayed_coord_transformer<T, V>::reallocate_col(size_t first, size_t last, Iter out) {
    if (first>last || last>row_index.size()) { 
        throw std::runtime_error("invalid row start/end indices");
    }
    
    const V& holding=tmp.vec;
    auto rIt=row_index.begin()+first, end=row_index.begin()+last;
    while (rIt!=end) {
        (*out)=holding[*rIt];
        ++out;
        ++rIt;
    }
    return;
}

/* Implementing methods for the 'delayed_matrix' class */

template<typename T, class V, class base_mat>
delayed_matrix<T, V, base_mat>::delayed_matrix(const Rcpp::RObject& incoming) : original(incoming), seed_ptr(nullptr) {
    check_DelayedMatrix(incoming);

    // Trying to generate the seed, if it's a valid object in itself.
    if (only_delayed_coord_changes(incoming)) {
        seed_ptr=generate_seed(incoming);
    }
        
    // If the seed is still NULL, we switch to a chunked matrix format.
    if (seed_ptr.get()==NULL) { 
        seed_ptr=generate_unknown_seed(incoming);
        transformer=delayed_coord_transformer<T, V>(seed_ptr.get());
    } else {
        transformer=delayed_coord_transformer<T, V>(incoming, seed_ptr.get());
    }

    nrow=transformer.get_nrow();
    ncol=transformer.get_ncol();
    return;
}

template<typename T, class V, class base_mat> 
delayed_matrix<T, V, base_mat>::~delayed_matrix() {}

template<typename T, class V, class base_mat> 
delayed_matrix<T, V, base_mat>::delayed_matrix(const delayed_matrix<T, V, base_mat>& other) : original(other.original), 
        seed_ptr(other.seed_ptr->clone()), transformer(other.transformer) {}

template<typename T, class V, class base_mat> 
delayed_matrix<T, V, base_mat>& delayed_matrix<T, V, base_mat>::operator=(const delayed_matrix<T, V, base_mat>& other) {
    original=other.original;
    seed_ptr=other.seed_ptr->clone();
    transformer=other.transformer;
    return *this;
}

template<typename T, class V, class base_mat>
template<class Iter>
void delayed_matrix<T, V, base_mat>::get_col(size_t c, Iter out, size_t first, size_t last) {
    transformer.get_col(seed_ptr.get(), c, out, first, last);
    return;
}

template<typename T, class V, class base_mat>
template<class Iter>
void delayed_matrix<T, V, base_mat>::get_row(size_t r, Iter out, size_t first, size_t last) {
    transformer.get_row(seed_ptr.get(), r, out, first, last);
    return;
}

template<typename T, class V, class base_mat>
T delayed_matrix<T, V, base_mat>::get(size_t r, size_t c) {
    return transformer.get(seed_ptr.get(), r, c);
}

template<typename T, class V, class base_mat> 
Rcpp::RObject delayed_matrix<T, V, base_mat>::yield() const {
    return original;
}

template<typename T, class V, class base_mat>
matrix_type delayed_matrix<T, V, base_mat>::get_matrix_type() const { 
    return DELAYED;
}

}
