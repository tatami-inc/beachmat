#ifndef BEACHMAT_CHARACTER_MATRIX_H
#define BEACHMAT_CHARACTER_MATRIX_H

#include "Input_matrix.h"

namespace beachmat { 

/* Virtual base class for character matrices. */

class character_matrix {
public:    
    character_matrix();
    virtual ~character_matrix();
    
    virtual size_t get_nrow() const=0;
    virtual size_t get_ncol() const=0;
    
    void get_row(size_t, Rcpp::StringVector::iterator); 
    virtual void get_row(size_t, Rcpp::StringVector::iterator, size_t, size_t)=0;

    void get_col(size_t, Rcpp::StringVector::iterator);
    virtual void get_col(size_t, Rcpp::StringVector::iterator, size_t, size_t)=0;

    virtual Rcpp::String get(size_t, size_t)=0;

    Rcpp::StringVector::iterator get_const_col(size_t, Rcpp::StringVector::iterator);
    virtual Rcpp::StringVector::iterator get_const_col(size_t, Rcpp::StringVector::iterator, size_t, size_t);

    const_col_indexed_info<Rcpp::StringVector> get_const_col_indexed(size_t, Rcpp::StringVector::iterator);
    virtual const_col_indexed_info<Rcpp::StringVector> get_const_col_indexed(size_t, Rcpp::StringVector::iterator, size_t, size_t);

    virtual std::unique_ptr<character_matrix> clone() const=0;

    virtual Rcpp::RObject yield () const=0;
    virtual matrix_type get_matrix_type() const=0;

private:
    Rcpp::IntegerVector indices; // needed for get_const_col_indexed for non-sparse matrices.
};

/* Advanced character matrix template */

template<class M>
class general_character_matrix : public character_matrix {
public:    
    general_character_matrix(const Rcpp::RObject& incoming) : mat(incoming) {}
    ~general_character_matrix() {}
  
    size_t get_nrow() const { return mat.get_nrow(); }
    size_t get_ncol() const { return mat.get_ncol(); }
 
    void get_row(size_t r, Rcpp::StringVector::iterator out, size_t first, size_t last) { return mat.get_row(r, out, first, last); }
    void get_col(size_t c, Rcpp::StringVector::iterator out, size_t first, size_t last) { return mat.get_col(c, out, first, last); }

    Rcpp::String get(size_t r, size_t c) { return mat.get(r, c); }

    std::unique_ptr<character_matrix> clone() const { return std::unique_ptr<character_matrix>(new general_character_matrix(*this)); }

    Rcpp::RObject yield () const { return mat.yield(); }
    matrix_type get_matrix_type() const { return mat.get_matrix_type(); }
protected:
    M mat;
};

/* Simple character matrix */

using simple_character_precursor=general_character_matrix<simple_matrix<Rcpp::String, Rcpp::StringVector> >;
class simple_character_matrix : public simple_character_precursor {
public:
    simple_character_matrix(const Rcpp::RObject& incoming);
    ~simple_character_matrix();
    Rcpp::StringVector::iterator get_const_col(size_t, Rcpp::StringVector::iterator, size_t, size_t);
    std::unique_ptr<character_matrix> clone() const;
};

/* RLE character matrix */

using Rle_character_matrix=general_character_matrix<Rle_matrix<Rcpp::String, Rcpp::StringVector> >;

/* HDF5Matrix */

class HDF5_character_matrix : public character_matrix {
public:    
    HDF5_character_matrix(const Rcpp::RObject&);
    ~HDF5_character_matrix();

    size_t get_nrow() const;
    size_t get_ncol() const;
 
    void get_row(size_t, Rcpp::StringVector::iterator, size_t, size_t);
    void get_col(size_t, Rcpp::StringVector::iterator, size_t, size_t);

    Rcpp::String get(size_t, size_t);

    std::unique_ptr<character_matrix> clone() const;

    Rcpp::RObject yield () const;
    matrix_type get_matrix_type() const;
protected:
    HDF5_matrix<char, STRSXP> mat; 
    H5::DataType str_type;

    size_t bufsize;
    std::vector<char> row_buf, col_buf, one_buf;
};

/* DelayedMatrix */

typedef delayed_matrix<Rcpp::String, Rcpp::StringVector, character_matrix> delayed_character_helper;

using delayed_character_matrix=general_character_matrix<delayed_character_helper>;

/* Unknown matrix type */

using unknown_character_matrix=general_character_matrix<unknown_matrix<Rcpp::String, Rcpp::StringVector> >;

/* Dispatcher */

std::unique_ptr<character_matrix> create_character_matrix(const Rcpp::RObject&);

}

/* Collected output definitions, so people only have to do #include "character_matrix.h" */

#include "character_output.h"

#endif
