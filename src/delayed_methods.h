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
delayed_coord_transformer<T, V>::delayed_coord_transformer(const Rcpp::List& net_subset, const Rcpp::LogicalVector& net_trans, M mat) : 
        transposed(false), byrow(false), bycol(false), delayed_nrow(mat->get_nrow()), delayed_ncol(mat->get_ncol()), 
        tmp(std::max(delayed_nrow, delayed_ncol)) {
   
    const size_t original_nrow(mat->get_nrow()), original_ncol(mat->get_ncol()); 

    if (net_subset.size()!=2) {
        throw std::runtime_error("subsetting list should be of length 2");
    }

    // Checking indices for rows.
    Rcpp::RObject rowdex(net_subset[0]);
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
    Rcpp::RObject coldex(net_subset[1]);
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

    // Checking transposition.
    if (net_trans.size()!=1) {
        throw std::runtime_error("transposition specifier should be of length 1");
    }
    transposed=net_trans[0];
    if (transposed) { // As the row/column indices refer to the matrix BEFORE transposition.
        std::swap(delayed_nrow, delayed_ncol);
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
    if (get_class(incoming)!="DelayedMatrix" || !incoming.isS4()) {
        throw std::runtime_error("input matrix should be a DelayedMatrix");
    }

    // Parsing the delayed operation structure.
    const Rcpp::Environment beachenv=Rcpp::Environment::namespace_env("beachmat");
    Rcpp::Function parser(beachenv["parseDelayedOps"]);
    Rcpp::List parse_out=parser(incoming);

    // Figuring if we can make use of the net subsetting/transposition state.
    if (parse_out.size()!=3) {
        throw std::runtime_error("output of beachmat:::parseDelayedOps should be a list of length 3");
    }
    seed_ptr=generate_seed(parse_out[2]);

    if (seed_ptr->get_matrix_type()!=UNKNOWN) {
        transformer=delayed_coord_transformer<T, V>(parse_out[0], parse_out[1], seed_ptr.get());
    } else {
        transformer=delayed_coord_transformer<T, V>(seed_ptr.get());
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
