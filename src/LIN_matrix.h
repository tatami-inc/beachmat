#ifndef BEACHMAT_LIN_MATRIX_H
#define BEACHMAT_LIN_MATRIX_H

#include "all_readers.h"

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

    /* We can't add a LogicalVector::iterator method because IntegerVector::iterator==LogicalVector::iterator
     * under the hood in Rcpp. The compiler then complains that overloading is not possible. Thus, for all 
     * references here to LogicalVector, we will consider the use of IntegerVector in its place.
     */

    void get_row(size_t, Rcpp::IntegerVector::iterator);
    void get_row(size_t, Rcpp::NumericVector::iterator);

    virtual void get_row(size_t, Rcpp::IntegerVector::iterator, size_t, size_t)=0;
    virtual void get_row(size_t, Rcpp::NumericVector::iterator, size_t, size_t)=0;

    void get_col(size_t, Rcpp::IntegerVector::iterator);
    void get_col(size_t, Rcpp::NumericVector::iterator);

    virtual void get_col(size_t, Rcpp::IntegerVector::iterator, size_t, size_t)=0;
    virtual void get_col(size_t, Rcpp::NumericVector::iterator, size_t, size_t)=0;

    virtual T get(size_t, size_t)=0;

    // Multi-row/column getters.
    void get_rows(Rcpp::IntegerVector::iterator, size_t, Rcpp::IntegerVector::iterator);
    void get_rows(Rcpp::IntegerVector::iterator ,size_t, Rcpp::NumericVector::iterator);

    virtual void get_rows(Rcpp::IntegerVector::iterator, size_t, Rcpp::IntegerVector::iterator, size_t, size_t)=0;
    virtual void get_rows(Rcpp::IntegerVector::iterator, size_t, Rcpp::NumericVector::iterator, size_t, size_t)=0;

    void get_cols(Rcpp::IntegerVector::iterator, size_t, Rcpp::IntegerVector::iterator);
    void get_cols(Rcpp::IntegerVector::iterator, size_t, Rcpp::NumericVector::iterator);

    virtual void get_cols(Rcpp::IntegerVector::iterator, size_t, Rcpp::IntegerVector::iterator, size_t, size_t)=0;
    virtual void get_cols(Rcpp::IntegerVector::iterator, size_t, Rcpp::NumericVector::iterator, size_t, size_t)=0;

    // Specialist getters.
    typename V::iterator get_const_col(size_t, typename V::iterator);
    virtual typename V::iterator get_const_col(size_t, typename V::iterator, size_t, size_t);

    const_col_indexed_info<V> get_const_col_indexed(size_t, typename V::iterator);
    virtual const_col_indexed_info<V> get_const_col_indexed(size_t, typename V::iterator, size_t, size_t);

    // Other methods.
    virtual std::unique_ptr<lin_matrix<T, V> > clone() const=0;

    virtual Rcpp::RObject yield() const=0;
    virtual matrix_type get_matrix_type() const=0;

private:
    Rcpp::IntegerVector indices; // needed for get_const_col_indexed for non-sparse matrices.
};

/* A general flavour for a LIN matrix */

template <typename T, class V, class RDR>
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

    void get_rows(Rcpp::IntegerVector::iterator, size_t, Rcpp::IntegerVector::iterator, size_t, size_t);
    void get_rows(Rcpp::IntegerVector::iterator, size_t, Rcpp::NumericVector::iterator, size_t, size_t);

    void get_cols(Rcpp::IntegerVector::iterator, size_t, Rcpp::IntegerVector::iterator, size_t, size_t);
    void get_cols(Rcpp::IntegerVector::iterator, size_t, Rcpp::NumericVector::iterator, size_t, size_t);

    std::unique_ptr<lin_matrix<T, V> > clone() const;

    Rcpp::RObject yield() const;
    matrix_type get_matrix_type() const;
protected:
    RDR reader;
};

/* Realization of the general flavour for a simple matrix */

template <typename T, class V>
using simple_lin_precursor=general_lin_matrix<T, V, simple_reader<T, V> >;

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
using dense_lin_precursor=general_lin_matrix<T, V, dense_reader<T, V> >;

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
using Csparse_lin_precursor=general_lin_matrix<T, V, Csparse_reader<T, V> >;

template <typename T, class V>
class Csparse_lin_matrix : public Csparse_lin_precursor<T, V> {
public:
    Csparse_lin_matrix(const Rcpp::RObject&);
    ~Csparse_lin_matrix();

    const_col_indexed_info<V> get_const_col_indexed(size_t, typename V::iterator, size_t, size_t);

    std::unique_ptr<lin_matrix<T, V> > clone() const;
};

/* HDF5Matrix of LINs */

template<typename T, int RTYPE>
class HDF5_lin_reader : public HDF5_reader<T, RTYPE> {
public:
    HDF5_lin_reader(const Rcpp::RObject&);
    ~HDF5_lin_reader();
    
    T get(size_t, size_t);
    
    void get_row(size_t, Rcpp::NumericVector::iterator, size_t, size_t);
    void get_row(size_t, Rcpp::IntegerVector::iterator, size_t, size_t);
    
    void get_col(size_t, Rcpp::NumericVector::iterator, size_t, size_t);
    void get_col(size_t, Rcpp::IntegerVector::iterator, size_t, size_t);

    void get_rows(Rcpp::IntegerVector::iterator, size_t, Rcpp::IntegerVector::iterator, size_t, size_t);
    void get_rows(Rcpp::IntegerVector::iterator, size_t, Rcpp::NumericVector::iterator, size_t, size_t);

    void get_cols(Rcpp::IntegerVector::iterator, size_t, Rcpp::IntegerVector::iterator, size_t, size_t);
    void get_cols(Rcpp::IntegerVector::iterator, size_t, Rcpp::NumericVector::iterator, size_t, size_t);
};

template <typename T, class V, int RTYPE>
using HDF5_lin_matrix=general_lin_matrix<T, V, HDF5_lin_reader<T, RTYPE> >;

/* DelayedMatrix of LINs */

template <typename T, class V>
using delayed_lin_reader=delayed_matrix<T, V, lin_matrix<T, V> >;

template <typename T, class V>
using delayed_lin_matrix=general_lin_matrix<T, V, delayed_lin_reader<T, V> >;

/* Unknown matrix of LINs */

template <typename T, class V>
using unknown_lin_matrix=general_lin_matrix<T, V, unknown_reader<T, V> >;

}

#include "LIN_methods_read.h"

#endif

