#ifndef BEACHMAT_LIN_MATRIX_H
#define BEACHMAT_LIN_MATRIX_H

#include "Rcpp.h"

#include "simple_reader.h"
#include "dense_reader.h"
#include "Csparse_reader.h"
#include "delayed_reader.h"
#include "unknown_reader.h"
#include "external_reader.h"
#include "../utils/utils.h"

#include <memory>

namespace beachmat { 

/***************************************************************** 
 * Virtual base class for LIN (logical/integer/numeric) matrices. 
 *****************************************************************/

template<typename T, class V>
class lin_matrix {
public:
    lin_matrix() = default;
    virtual ~lin_matrix() = default;
    lin_matrix(const lin_matrix&) = default;
    lin_matrix& operator=(const lin_matrix&) = default;
    lin_matrix(lin_matrix&&) = default;
    lin_matrix& operator=(lin_matrix&&) = default;
    
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

    virtual std::string get_class() const=0;

    virtual std::string get_package() const=0;

private:
    Rcpp::IntegerVector indices; // needed for get_const_col_indexed for non-sparse matrices.
};

/* A general flavour for a LIN matrix */

template <typename T, class V, class RDR>
class general_lin_matrix : public lin_matrix<T, V> {
public:    
    general_lin_matrix(const Rcpp::RObject&);
    ~general_lin_matrix() = default;
    general_lin_matrix(const general_lin_matrix&) = default;
    general_lin_matrix& operator=(const general_lin_matrix&) = default;
    general_lin_matrix(general_lin_matrix&&) = default;
    general_lin_matrix& operator=(general_lin_matrix&&) = default;
    
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

    std::string get_class() const { return reader.get_class(); }

    std::string get_package() const { return reader.get_package(); }
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
    ~simple_lin_matrix() = default;
    simple_lin_matrix(const simple_lin_matrix&) = default;
    simple_lin_matrix& operator=(const simple_lin_matrix&) = default;
    simple_lin_matrix(simple_lin_matrix&&) = default;
    simple_lin_matrix& operator=(simple_lin_matrix&&) = default;
    
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
    ~dense_lin_matrix() = default;
    dense_lin_matrix(const dense_lin_matrix&) = default;
    dense_lin_matrix& operator=(const dense_lin_matrix&) = default;
    dense_lin_matrix(dense_lin_matrix&&) = default;
    dense_lin_matrix& operator=(dense_lin_matrix&&) = default;

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
    ~Csparse_lin_matrix() = default;
    Csparse_lin_matrix(const Csparse_lin_matrix&) = default;
    Csparse_lin_matrix& operator=(const Csparse_lin_matrix&) = default;
    Csparse_lin_matrix(Csparse_lin_matrix&&) = default;
    Csparse_lin_matrix& operator=(Csparse_lin_matrix&&) = default;

    const_col_indexed_info<V> get_const_col_indexed(size_t, typename V::iterator, size_t, size_t);

    std::unique_ptr<lin_matrix<T, V> > clone() const;
};

/* DelayedMatrix of LINs */

template <typename T, class V>
using delayed_lin_reader=delayed_reader<T, V, lin_matrix<T, V> >;

template <typename T, class V>
using delayed_lin_matrix=general_lin_matrix<T, V, delayed_lin_reader<T, V> >;

/* Unknown matrix of LINs */

template <typename T, class V>
using unknown_lin_matrix=general_lin_matrix<T, V, unknown_reader<T, V> >;

/* External matrix of LINs */

template <typename T, class V>
using external_lin_precursor=general_lin_matrix<T, V, external_lin_reader<T, V> >;

template <typename T, class V>
class external_lin_matrix : public external_lin_precursor<T, V> {
public:
    external_lin_matrix(const Rcpp::RObject&);
    ~external_lin_matrix() = default;
    external_lin_matrix(const external_lin_matrix&) = default;
    external_lin_matrix& operator=(const external_lin_matrix&) = default;
    external_lin_matrix(external_lin_matrix&&) = default;
    external_lin_matrix& operator=(external_lin_matrix&&) = default;
    
    typename V::iterator get_const_col(size_t, typename V::iterator, size_t, size_t);

    const_col_indexed_info<V> get_const_col_indexed(size_t, typename V::iterator, size_t, size_t);

    std::unique_ptr<lin_matrix<T, V> > clone() const;
};

}

#include "LIN_methods.h"

#endif

