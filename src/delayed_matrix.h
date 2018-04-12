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
    delayed_coord_transformer();

    template<class M>
    delayed_coord_transformer(M);

    template<class M>
    delayed_coord_transformer(const Rcpp::List&, const Rcpp::LogicalVector&, M);

    template<class M, class Iter>
    void get_row(M, size_t, Iter, size_t, size_t);
    
    template<class M, class Iter>
    void get_col(M, size_t, Iter, size_t, size_t);

    template<class M>
    T get(M, size_t, size_t);

    size_t get_nrow() const;
    size_t get_ncol() const;
private:
    std::vector<int> row_index, col_index;
    bool transposed, byrow, bycol;
    size_t delayed_nrow, delayed_ncol;

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

/* The 'delayed_matrix' class, which wraps the coord_transformer class. */

template<typename T, class V, class base_mat>
class delayed_matrix : public any_matrix { 
public:
    delayed_matrix(const Rcpp::RObject&);
    ~delayed_matrix();
    delayed_matrix(const delayed_matrix&);
    delayed_matrix& operator=(const delayed_matrix&);
    delayed_matrix(delayed_matrix&&) = default;             // Move constructor
    delayed_matrix& operator=(delayed_matrix&&) = default;  // Move assignment constructor
    
    T get(size_t, size_t);

    template <class Iter>
    void get_row(size_t, Iter, size_t, size_t);

    template <class Iter>
    void get_col(size_t, Iter, size_t, size_t);

    Rcpp::RObject yield() const;
    matrix_type get_matrix_type() const;
private:
    Rcpp::RObject original;
    std::unique_ptr<base_mat> seed_ptr;
    delayed_coord_transformer<T, V> transformer;

    // Specialized function for each realized matrix type.
    static std::unique_ptr<base_mat> generate_seed(Rcpp::RObject);
};

}

#include "delayed_methods.h"

#endif 
