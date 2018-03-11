#ifndef BEACHMAT_LIN_MATRIX_H
#define BEACHMAT_LIN_MATRIX_H

#include "Input_matrix.h"

namespace beachmat { 

/***************************************************************** 
 * Virtual base class for LIN (logical/integer/numeric) matrices. 
 *****************************************************************/

template<typename T, class V>
class lin_matrix {
public:
    lin_matrix();
    virtual ~lin_matrix();
    
    virtual size_t get_nrow() const=0;
    virtual size_t get_ncol() const=0;

    void get_row(size_t, Rcpp::IntegerVector::iterator);
    void get_row(size_t, Rcpp::NumericVector::iterator);

    /* We can't add a LogicalVector::iterator method because IntegerVector::iterator==LogicalVector::iterator
     * under the hood in Rcpp. The compiler then complains that overloading is not possible. Thus, for all 
     * references here to LogicalVector, we will consider the use of IntegerVector in its place.
     */

    virtual void get_row(size_t, Rcpp::IntegerVector::iterator, size_t, size_t)=0;
    virtual void get_row(size_t, Rcpp::NumericVector::iterator, size_t, size_t)=0;

    void get_col(size_t, Rcpp::IntegerVector::iterator);
    void get_col(size_t, Rcpp::NumericVector::iterator);

    virtual void get_col(size_t, Rcpp::IntegerVector::iterator, size_t, size_t)=0;
    virtual void get_col(size_t, Rcpp::NumericVector::iterator, size_t, size_t)=0;

    virtual T get(size_t, size_t)=0;

    typename V::iterator get_const_col(size_t, typename V::iterator);
    virtual typename V::iterator get_const_col(size_t, typename V::iterator, size_t, size_t);

    typedef std::tuple<size_t, Rcpp::IntegerVector::iterator, typename V::iterator> const_col_indexed_info;
    const_col_indexed_info get_const_col_indexed(size_t, typename V::iterator);
    virtual const_col_indexed_info get_const_col_indexed(size_t, typename V::iterator, size_t, size_t);

    virtual std::unique_ptr<lin_matrix<T, V> > clone() const=0;

    virtual Rcpp::RObject yield() const=0;
    virtual matrix_type get_matrix_type() const=0;

private:
    Rcpp::IntegerVector indices; // needed for get_const_col_indexed for non-sparse matrices.
};

/* A general flavour for a LIN matrix */

template <typename T, class V, class M>
class general_lin_matrix : public lin_matrix<T, V> {
public:    
    general_lin_matrix(const Rcpp::RObject&);
    ~general_lin_matrix();
    
    size_t get_nrow() const;
    size_t get_ncol() const;

    void get_col(size_t,  Rcpp::IntegerVector::iterator, size_t, size_t);
    void get_col(size_t,  Rcpp::NumericVector::iterator, size_t, size_t);

    void get_row(size_t,  Rcpp::IntegerVector::iterator, size_t, size_t);
    void get_row(size_t,  Rcpp::NumericVector::iterator, size_t, size_t);

    T get(size_t, size_t);

    std::unique_ptr<lin_matrix<T, V> > clone() const;

    Rcpp::RObject yield() const;
    matrix_type get_matrix_type() const;
protected:
    M mat;
};

/* Realization of the general flavour for a simple matrix */

template <typename T, class V>
using simple_lin_precursor=general_lin_matrix<T, V, simple_matrix<T, V> >;

template <typename T, class V>
class simple_lin_matrix : public simple_lin_precursor<T, V> {
public:
    simple_lin_matrix(const Rcpp::RObject&);
    ~simple_lin_matrix();
    
    typename V::iterator get_const_col(size_t, typename V::iterator, size_t, size_t);

    std::unique_ptr<lin_matrix<T, V> > clone() const;
};

/* Realization of the general flavour for a dense matrix */

template <typename T, class V>
using dense_lin_precursor=general_lin_matrix<T, V, dense_matrix<T, V> >;

template <typename T, class V>
class dense_lin_matrix : public dense_lin_precursor<T, V> {
public:
    dense_lin_matrix(const Rcpp::RObject&);
    ~dense_lin_matrix();
    
    typename V::iterator get_const_col(size_t, typename V::iterator, size_t, size_t);

    std::unique_ptr<lin_matrix<T, V> > clone() const;
};

/* Realization of the general flavour for a C-sparse matrix */

template <typename T, class V>
using Csparse_lin_precursor=general_lin_matrix<T, V, Csparse_matrix<T, V> >;

template <typename T, class V>
class Csparse_lin_matrix : public Csparse_lin_precursor<T, V> {
public:
    Csparse_lin_matrix(const Rcpp::RObject&);
    ~Csparse_lin_matrix();

    typename lin_matrix<T, V>::const_col_indexed_info get_const_col_indexed(size_t, typename V::iterator, size_t, size_t);

    std::unique_ptr<lin_matrix<T, V> > clone() const;
};

/* Realization of the general flavour for other matrices */

template <typename T, class V>
using Psymm_lin_matrix=general_lin_matrix<T, V, Psymm_matrix<T, V> >;

template <typename T, class V>
using Rle_lin_matrix=general_lin_matrix<T, V, Rle_matrix<T, V> >;

/* HDF5Matrix of LINs */

template<typename T, class V, int RTYPE>
class HDF5_lin_matrix : public lin_matrix<T, V> {
public:
    HDF5_lin_matrix(const Rcpp::RObject&);
    ~HDF5_lin_matrix();

    size_t get_nrow() const;
    size_t get_ncol() const;

    void get_col(size_t, Rcpp::IntegerVector::iterator, size_t, size_t);
    void get_col(size_t, Rcpp::NumericVector::iterator, size_t, size_t);

    void get_row(size_t, Rcpp::IntegerVector::iterator, size_t, size_t);
    void get_row(size_t, Rcpp::NumericVector::iterator, size_t, size_t);

    T get(size_t, size_t);

    std::unique_ptr<lin_matrix<T, V> > clone() const;

    Rcpp::RObject yield() const;
    matrix_type get_matrix_type() const;
protected:
    HDF5_matrix<T, RTYPE> mat;
};

/* DelayedMatrix of LINs */

template<typename T, class V>
class delayed_lin_matrix : public lin_matrix<T, V> {
public:
    delayed_lin_matrix(const Rcpp::RObject&);
    ~delayed_lin_matrix();
    delayed_lin_matrix(const delayed_lin_matrix&);
    delayed_lin_matrix& operator=(const delayed_lin_matrix&);
    delayed_lin_matrix(delayed_lin_matrix&&) = default;             // Move constructor
    delayed_lin_matrix& operator=(delayed_lin_matrix&&) = default;  // Move assignment constructor
    
    size_t get_nrow() const;
    size_t get_ncol() const;
    
    void get_col(size_t, Rcpp::IntegerVector::iterator, size_t, size_t);
    void get_col(size_t, Rcpp::NumericVector::iterator, size_t, size_t);

    void get_row(size_t, Rcpp::IntegerVector::iterator, size_t, size_t);
    void get_row(size_t, Rcpp::NumericVector::iterator, size_t, size_t);

    T get(size_t, size_t);

    std::unique_ptr<lin_matrix<T, V> > clone() const;

    Rcpp::RObject yield() const;
    matrix_type get_matrix_type() const;
   
private:
    Rcpp::RObject original;
    std::unique_ptr<lin_matrix<T, V> > seed_ptr;
    delayed_coord_transformer<T, V> transformer;
    static std::unique_ptr<lin_matrix<T, V> > generate_seed(Rcpp::RObject);

    /* This class allows chunked extraction from a 'delayed_matrix' instance,
     * mimicking a unique pointer to other top-level classes.
     */
    using enslaved_precursor=general_lin_matrix<T, V, delayed_matrix<T, V> >;
    class enslaved : public enslaved_precursor {
    public:
        enslaved(const Rcpp::RObject&);
        ~enslaved();
//        typename V::iterator get_const_col(size_t, typename V::iterator, size_t, size_t);
        std::unique_ptr<lin_matrix<T, V> > clone() const;
    };
};

}

#include "LIN_methods.h"

#endif

