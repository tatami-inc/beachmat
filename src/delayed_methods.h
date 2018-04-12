#include "delayed_matrix.h"
#include "unknown_matrix.h"

namespace beachmat {

/* Implementing methods for the 'delayed_coord_transformer' class */

template<typename T, class V>
delayed_coord_transformer<T, V>::delayed_coord_transformer() {}

template<typename T, class V>
template<class M>
delayed_coord_transformer<T, V>::delayed_coord_transformer(M mat) : delayed_nrow(mat->get_nrow()), delayed_ncol(mat->get_ncol()) {}

template<typename T, class V>
template<class M>
delayed_coord_transformer<T, V>::delayed_coord_transformer(const Rcpp::List& net_subset, const Rcpp::LogicalVector& net_trans, M mat) :
        delayed_nrow(mat->get_nrow()), delayed_ncol(mat->get_ncol()), tmp(std::max(delayed_nrow, delayed_ncol)) {
   
    const size_t original_nrow(mat->get_nrow()), original_ncol(mat->get_ncol()); 

    if (net_subset.size()!=2) {
        throw std::runtime_error("subsetting list should be of length 2");
    }

    // Checking indices for rows.
    obtain_indices(net_subset[0], original_nrow, byrow, delayed_nrow, row_index);

    // Checking indices for columns.
    obtain_indices(net_subset[1], original_ncol, bycol, delayed_ncol, col_index);

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
void delayed_coord_transformer<T, V>::obtain_indices(const Rcpp::RObject& subset_in, size_t original_dim,
        bool& affected, size_t& delayed_dim, std::vector<size_t>& subset_out) {
    // This function simply converts the row or column subset indices to 0-indexed form,
    // setting affected=false if there are no subset indices or if the subset indices are 1:original_dim.
    // Note that delayed_dim is also reset to the length fo the subset index vector.

    affected=!subset_in.isNULL();
    if (!affected){ 
        return;
    }

    // Coercing the subset indices to zero-indexed size_t's.
    if (subset_in.sexp_type()!=INTSXP) {
        throw std::runtime_error("index vector should be integer");
    }

    Rcpp::IntegerVector idx(subset_in);
    delayed_dim=idx.size();
    subset_out.reserve(delayed_dim);

    for (const auto& i : idx) {
        if (i < 1 || i > original_dim) {
            throw std::runtime_error("delayed subset indices are out of range");
        }
        subset_out.push_back(i-1);
    }

    // If the indices are all consecutive from 0 to N-1, we turn off 'affected'. 
    if (delayed_dim 
            && delayed_dim==original_dim
            && subset_out.front()==0 
            && subset_out.back()+1==delayed_dim) {

        size_t count=0;
        affected=false;
        for (auto i : subset_out) {
            if (i!=count) {
                affected=true;
                break;
            }
            ++count;
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
            reallocate_col(mat, r, first, last, out);
        } else {
            mat->get_col(r, out, first, last);
        }
    } else {
        r=transform_row(r);

        // Row extraction, first/last refer to columns.
        if (bycol) {
            reallocate_row(mat, r, first, last, out);
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
            reallocate_row(mat, c, first, last, out);
        } else {
            mat->get_row(c, out, first, last);
        }
    } else {
        c=transform_col(c);

        // Column extraction, first/last refer to rows.
        if (byrow) {
            reallocate_col(mat, c, first, last, out);
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
void delayed_coord_transformer<T, V>::prepare_reallocation(size_t first, size_t last, 
        size_t& old_first, size_t& old_last, size_t& min_index, size_t& max_index, 
        const std::vector<size_t>& indices, const std::string& msg) {

    if (first>last || last>indices.size()) { 
        throw_custom_error("invalid ", msg, " start/end indices");
    }

    if (old_first!=first || old_last!=last) {
        old_first=first;
        old_last=last;
        if (first!=last) {
            min_index=*std::min_element(indices.begin()+first, indices.begin()+last);
            max_index=*std::max_element(indices.begin()+first, indices.begin()+last)+1;
        } else {
            // Avoid problems with max/min of zero-length vectors.
            min_index=0;
            max_index=0;
        }
    }

    return;
}

template<typename T, class V>
template<class M, class Iter>
void delayed_coord_transformer<T, V>::reallocate_row(M mat, size_t r, size_t first, size_t last, Iter out) {
    // Yes, the use of *_col_* variables for row reallocation is intentional!
    // This is because we're talking about the columns in the extracted row.
    prepare_reallocation(first, last, old_col_first, old_col_last, 
            min_col_index, max_col_index, col_index, "column");

    V& holding=tmp.vec;
    mat->get_row(r, holding.begin(), min_col_index, max_col_index);
    auto cIt=col_index.begin()+first, end=col_index.begin()+last;
    while (cIt!=end) {
        (*out)=holding[*cIt - min_col_index];
        ++out;
        ++cIt;
    }
    return;
}

template<typename T, class V>
template<class M, class Iter>
void delayed_coord_transformer<T, V>::reallocate_col(M mat, size_t c, size_t first, size_t last, Iter out) {
    // Yes, the use of *_row_* variables for column reallocation is intentional!
    // This is because we're talking about the rows in the extracted column.
    prepare_reallocation(first, last, old_row_first, old_row_last, 
            min_row_index, max_row_index, row_index, "row");

    V& holding=tmp.vec;
    mat->get_col(c, holding.begin(), min_row_index, max_row_index);
    auto rIt=row_index.begin()+first, end=row_index.begin()+last;
    while (rIt!=end) {
        (*out)=holding[*rIt - min_row_index];
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
