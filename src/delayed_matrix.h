#ifndef BEACHMAT_DELAYED_MATRIX_H
#define BEACHMAT_DELAYED_MATRIX_H
#include "beachmat.h"
#include "utils.h"

namespace beachmat {

/* The 'delayed_coord_transformer' class is intended to allow direct data access
 * from the underlying seed matrix in a DelayedMatrix class, while transforming
 * the extracted values to account for row/column subsetting or transposition.
 * This avoids the need to realize any subset of the matrix, as would be necessary
 * for more general delayed operations (see the 'delayed_matrix' class below).
 */

template<typename T, class V>
class delayed_coord_transformer {
public:    
    delayed_coord_transformer(const Rcpp::RObject&);

    template<class M>
    void set_dim(M mat);

    template<class M, class Iter>
    void get_row(M, size_t, Iter, size_t, size_t);
    
    template<class M, class Iter>
    void get_col(M, size_t, Iter, size_t, size_t);

    template<class M>
    T get(M, size_t, size_t);

    size_t get_nrow() const;
    size_t get_ncol() const;
    bool has_unmodified_values() const;
private:
    std::vector<int> row_index, col_index;
    bool transposed, byrow, bycol;
    bool coord_changes_only;

    size_t original_nrow, original_ncol, delayed_nrow, delayed_ncol;

    /* Making a copyable vector to save myself having to write copy constructors for the entire transformer.
     * This is necessary as we need an internal Rcpp::Vector to extract from the underlying seed
     * prior to subsetting, but such vectors are not actually default-copied between instances.
     */
    struct copyable_holder {
        copyable_holder(size_t n=0) : vec(n) {}
        ~copyable_holder() {};
        copyable_holder(const copyable_holder& other) : vec(Rcpp::clone(other.vec)) {}
        copyable_holder& operator=(const copyable_holder& other) { 
            vec=Rcpp::clone(other.vec); 
            return *this;
        }
        copyable_holder(copyable_holder&&) = default;
        copyable_holder& operator=(copyable_holder&&) = default;
        V vec;
    };
    copyable_holder tmp;

    size_t transform_row(size_t) const;
    size_t transform_col(size_t) const;

    template<class Iter>
    void reallocate_row(size_t, size_t, Iter out);
    template<class Iter>
    void reallocate_col(size_t, size_t, Iter out);
};

/* The 'delayed_matrix' class will realize chunks of the input DelayedMatrix
 * upon request from any calling function. This is necessary to avoid 
 * reimplementing arbitrary delayed R operations in C++.
 */

template<typename T, class V>
class delayed_matrix : public any_matrix {
public:    
    delayed_matrix(const Rcpp::RObject&);
    ~delayed_matrix();

    T get(size_t, size_t);

    template <class Iter> 
    void get_row(size_t, Iter, size_t, size_t);

    template <class Iter>
    void get_col(size_t, Iter, size_t, size_t);

    typename V::iterator get_const_col(size_t, size_t, size_t);
    
    Rcpp::RObject yield() const;
    matrix_type get_matrix_type () const;
private:
    Rcpp::RObject original;

    Rcpp::Environment beachenv;
    Rcpp::Function realizer_row, realizer_col;
    
    V storage;
    void update_storage_by_row(size_t);
    void update_storage_by_col(size_t);

    Rcpp::IntegerVector row_indices, col_indices;
    int chunk_nrow, chunk_ncol;
};

}

#include "delayed_methods.h"

#endif 
