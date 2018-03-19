#ifndef BEACHMAT_LIN_OUTPUT_H
#define BEACHMAT_LIN_OUTPUT_H

#include "Output_matrix.h"
#include "output_param.h"

namespace beachmat { 

/************************************************************************
 * Virtual base class for LIN (logical/integer/numeric) output matrices.
 ************************************************************************/

template<typename T, class V>
class lin_output {
public:
    lin_output();
    virtual ~lin_output();
    
    // Getters:
    virtual size_t get_nrow() const=0;
    virtual size_t get_ncol() const=0;

    void get_row(size_t, Rcpp::IntegerVector::iterator);
    void get_row(size_t, Rcpp::NumericVector::iterator);

    virtual void get_row(size_t, Rcpp::IntegerVector::iterator, size_t, size_t)=0;
    virtual void get_row(size_t, Rcpp::NumericVector::iterator, size_t, size_t)=0;

    void get_col(size_t, Rcpp::IntegerVector::iterator);
    void get_col(size_t, Rcpp::NumericVector::iterator);

    virtual void get_col(size_t, Rcpp::IntegerVector::iterator, size_t, size_t)=0;
    virtual void get_col(size_t, Rcpp::NumericVector::iterator, size_t, size_t)=0;

    typename V::iterator get_const_col(size_t, typename V::iterator);
    virtual typename V::iterator get_const_col(size_t, typename V::iterator, size_t, size_t);

    virtual T get(size_t, size_t)=0;

    // Setters:
    void set_row(size_t, Rcpp::IntegerVector::iterator);
    void set_row(size_t, Rcpp::NumericVector::iterator);

    virtual void set_row(size_t, Rcpp::IntegerVector::iterator, size_t, size_t)=0;
    virtual void set_row(size_t, Rcpp::NumericVector::iterator, size_t, size_t)=0;

    void set_col(size_t, Rcpp::IntegerVector::iterator);
    void set_col(size_t, Rcpp::NumericVector::iterator);

    virtual void set_col(size_t, Rcpp::IntegerVector::iterator, size_t, size_t)=0;
    virtual void set_col(size_t, Rcpp::NumericVector::iterator, size_t, size_t)=0;

    virtual void set_col_indexed(size_t, const const_col_indexed_info<Rcpp::IntegerVector>&)=0;
    virtual void set_col_indexed(size_t, const const_col_indexed_info<Rcpp::NumericVector>&)=0;

    virtual void set(size_t, size_t, T)=0;

    // Other methods:
    virtual Rcpp::RObject yield()=0;

    virtual std::unique_ptr<lin_output<T, V> > clone() const=0;

    virtual matrix_type get_matrix_type() const=0;
private:
    Rcpp::IntegerVector indices; // needed for get_const_col_indexed.
};

/* General output */

template<typename T, class V, class M>
class general_lin_output : public lin_output<T, V> {
public:
    general_lin_output(size_t, size_t);
    ~general_lin_output();

    size_t get_nrow() const;
    size_t get_ncol() const;

    void get_row(size_t, Rcpp::IntegerVector::iterator, size_t, size_t);
    void get_row(size_t, Rcpp::NumericVector::iterator, size_t, size_t);

    void get_col(size_t, Rcpp::IntegerVector::iterator, size_t, size_t);
    void get_col(size_t, Rcpp::NumericVector::iterator, size_t, size_t);

    T get(size_t, size_t);

    void set_row(size_t, Rcpp::IntegerVector::iterator, size_t, size_t);
    void set_row(size_t, Rcpp::NumericVector::iterator, size_t, size_t);

    void set_col(size_t, Rcpp::IntegerVector::iterator, size_t, size_t);
    void set_col(size_t, Rcpp::NumericVector::iterator, size_t, size_t);

    void set_col_indexed(size_t, const const_col_indexed_info<Rcpp::IntegerVector>&);
    void set_col_indexed(size_t, const const_col_indexed_info<Rcpp::NumericVector>&);

    void set(size_t, size_t, T);

    Rcpp::RObject yield();

    std::unique_ptr<lin_output<T, V> > clone() const;

    matrix_type get_matrix_type() const;
protected:
    M mat;
};

/* Simple LIN output */

template<typename T, class V>
using simple_lin_output_precursor=general_lin_output<T, V, simple_output<T, V> >;

template<typename T, class V>
class simple_lin_output : public simple_lin_output_precursor<T, V> {
public:
    simple_lin_output(size_t, size_t);
    ~simple_lin_output();
    typename V::iterator get_const_col(size_t, typename V::iterator, size_t, size_t);
    std::unique_ptr<lin_output<T, V> > clone() const;
};

/* Sparse LIN output */

template<typename T, class V>
using sparse_lin_output_precursor=general_lin_output<T, V, Csparse_output<T, V> >;

template<typename T, class V>
class sparse_lin_output : public sparse_lin_output_precursor<T, V> {
public:
    sparse_lin_output(size_t, size_t);
    ~sparse_lin_output();
    std::unique_ptr<lin_output<T, V> > clone() const;
};

/* HDF5 LIN output */

template<typename T, class V>
class HDF5_lin_output : public lin_output<T, V> {
public:
    HDF5_lin_output(size_t, size_t, 
            size_t=output_param::DEFAULT_CHUNKDIM, 
            size_t=output_param::DEFAULT_CHUNKDIM, 
            int=output_param::DEFAULT_COMPRESS);
    ~HDF5_lin_output();

    size_t get_nrow() const;
    size_t get_ncol() const;

    void get_row(size_t, Rcpp::IntegerVector::iterator, size_t, size_t);
    void get_row(size_t, Rcpp::NumericVector::iterator, size_t, size_t);

    void get_col(size_t, Rcpp::IntegerVector::iterator, size_t, size_t);
    void get_col(size_t, Rcpp::NumericVector::iterator, size_t, size_t);

    T get(size_t, size_t);

    void set_row(size_t, Rcpp::IntegerVector::iterator, size_t, size_t);
    void set_row(size_t, Rcpp::NumericVector::iterator, size_t, size_t);

    void set_col(size_t, Rcpp::IntegerVector::iterator, size_t, size_t);
    void set_col(size_t, Rcpp::NumericVector::iterator, size_t, size_t);

    void set(size_t, size_t, T);

    void set_col_indexed(size_t, const const_col_indexed_info<Rcpp::IntegerVector>&);
    void set_col_indexed(size_t, const const_col_indexed_info<Rcpp::NumericVector>&);

    Rcpp::RObject yield();

    std::unique_ptr<lin_output<T, V> > clone() const;

    matrix_type get_matrix_type() const;
protected:
    HDF5_output<T, V> mat;
};

}

#include "LIN_outfun.h"

#endif

