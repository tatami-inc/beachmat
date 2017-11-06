#ifndef BEACHMAT_DELAYED_MATRIX_H
#define BEACHMAT_DELAYED_MATRIX_H
#include "beachmat.h"
#include "utils.h"

namespace beachmat {

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

    // Making a copyable vector to save myself having to write copy constructors for the transformer.
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

}

#include "delayed_methods.h"

#endif 
