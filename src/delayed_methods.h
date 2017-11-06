#include "delayed_matrix.h"

namespace beachmat {

template<typename T, class V>
delayed_coord_transformer<T, V>::delayed_coord_transformer(const Rcpp::RObject& in) : original_nrow(0), original_ncol(0), delayed_nrow(0), delayed_ncol(0) {
    if (get_class(in)!="DelayedMatrix") { 
        throw std::runtime_error("input matrix should be a DelayedMatrix");
    }

    // Checking if there's only coordinate changes.
    Rcpp::List delayed_ops(get_safe_slot(in, "delayed_ops"));
    coord_changes_only=(delayed_ops.size()==0);

    // Checking indices.
    Rcpp::List indices(get_safe_slot(in, "index"));

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
    }

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
    }
    
    // CHecking transposition.
    Rcpp::LogicalVector is_trans(get_safe_slot(in, "is_transposed"));
    if (is_trans.size()!=1) { 
        throw std::runtime_error("transposition specifier should be a logical scalar");
    }
    transposed=is_trans[0];
    return;
}

template<typename T, class V>
template<class M>
void delayed_coord_transformer<T, V>::set_dim(M mat) {
    delayed_nrow=original_nrow=mat->get_nrow();
    delayed_ncol=original_ncol=mat->get_ncol();
    tmp.vec=V(std::max(original_ncol, original_nrow));

    if (byrow) {
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

    if (bycol) { 
        delayed_ncol=col_index.size();

        // Same for bycol.
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

    if (transposed) {
        std::swap(delayed_nrow, delayed_ncol);
    }
    return;
};

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

template<typename T, class V>
bool delayed_coord_transformer<T, V>::has_unmodified_values() const {
    return coord_changes_only;
}

}
