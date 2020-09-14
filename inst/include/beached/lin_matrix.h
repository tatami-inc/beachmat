#ifndef BEACHMAT_LIN_MATRIX_H
#define BEACHMAT_LIN_MATRIX_H

#include "Rcpp.h"
#include "ordinary_reader.h"
#include "Csparse_reader.h"
#include "utils.h"

namespace beachmat {

/*** Base classes ***/

class lin_matrix {
public:
    lin_matrix() {}

    virtual ~lin_matrix() = default;
    lin_matrix(const lin_matrix&) = default;
    lin_matrix& operator=(const lin_matrix&) = default;
    lin_matrix(lin_matrix&&) = default;
    lin_matrix& operator=(lin_matrix&&) = default;

    virtual const int* get_col(size_t, int*, size_t, size_t) = 0;

    virtual const int* get_row(size_t, int*, size_t, size_t) = 0;

    virtual const double* get_col(size_t, double*, size_t, size_t) = 0;

    virtual const double* get_row(size_t, double*, size_t, size_t) = 0;

    const int* get_col(size_t r, int* work) {
        return get_col(r, work, 0, nrow);
    }

    const int* get_row(size_t c, int* work) {
        return get_row(c, work, 0, ncol);        
    }

    const double* get_col(size_t r, double* work) {
        return get_col(r, work, 0, nrow);
    }

    const double* get_row(size_t c, double* work) {
        return get_row(c, work, 0, ncol);        
    }

    size_t get_nrow() const { return nrow; }

    size_t get_ncol() const { return ncol; }
protected:
    size_t nrow=0, ncol=0;
};

class sparse_lin_matrix : public lin_matrix {
public:
    sparse_lin_matrix() {}

    ~sparse_lin_matrix() = default;
    sparse_lin_matrix(const sparse_lin_matrix&) = default;
    sparse_lin_matrix& operator=(const sparse_lin_matrix&) = default;
    sparse_lin_matrix(sparse_lin_matrix&&) = default;
    sparse_lin_matrix& operator=(sparse_lin_matrix&&) = default;

    virtual sparse_index<const int*, int> get_row(size_t, int*, int*, size_t, size_t) = 0;

    virtual sparse_index<const int*, int> get_col(size_t, int*, int*, size_t, size_t) = 0;

    virtual sparse_index<const double*, int> get_row(size_t, double*, int*, size_t, size_t) = 0;

    virtual sparse_index<const double*, int> get_col(size_t, double*, int*, size_t, size_t) = 0;

    sparse_index<const int*, int> get_col(size_t c, int* work_x, int* work_i) {
        return get_col(c, work_x, work_i, 0, this->nrow);
    }

    sparse_index<const int*, int> get_row(size_t r, int* work_x, int* work_i) {
        return get_row(r, work_x, work_i, 0, this->ncol);
    }

    sparse_index<const double*, int> get_col(size_t c, double* work_x, int* work_i) {
        return get_col(c, work_x, work_i, 0, this->nrow);
    }

    sparse_index<const double*, int> get_row(size_t r, double* work_x, int* work_i) {
        return get_row(r, work_x, work_i, 0, this->ncol);
    }
};

/*** Ordinary matrix ***/

template <class V>
class ordinary_matrix : public lin_matrix {
public:
    ordinary_matrix(Rcpp::RObject mat) : reader(mat) {
        this->nrow = reader.get_nrow();
        this->ncol = reader.get_ncol();
        return;
    }

    ~ordinary_matrix() = default;
    ordinary_matrix(const ordinary_matrix&) = default;
    ordinary_matrix& operator=(const ordinary_matrix&) = default;
    ordinary_matrix(ordinary_matrix&&) = default;
    ordinary_matrix& operator=(ordinary_matrix&&) = default;

    const int* get_col(size_t c, int* work, size_t first, size_t last);

    const int* get_row(size_t r, int* work, size_t first, size_t last) {
        reader.get_row(r, work, first, last);       
        return work;
    }

    const double* get_col(size_t c, double* work, size_t first, size_t last);

    const double* get_row(size_t r, double* work, size_t first, size_t last) {
        reader.get_row(r, work, first, last); 
        return work;
    }
private:
    ordinary_reader<V> reader;
};

using ordinary_integer_matrix = ordinary_matrix<Rcpp::IntegerVector>;

template <>
const int* ordinary_integer_matrix::get_col(size_t c, int* work, size_t first, size_t last) {
    return reader.get_col(c, first, last);
}

template <>
const double* ordinary_integer_matrix::get_col(size_t c, double* work, size_t first, size_t last) {
    auto out = reader.get_col(c, first, last);
    std::copy(out, out + last - first, work);
    return work;
}

using ordinary_logical_matrix = ordinary_matrix<Rcpp::LogicalVector>;

template <>
const int* ordinary_logical_matrix::get_col(size_t c, int* work, size_t first, size_t last) {
    return reader.get_col(c, first, last);
}

template <>
const double* ordinary_logical_matrix::get_col(size_t c, double* work, size_t first, size_t last) {
    auto out = reader.get_col(c, first, last);
    std::copy(out, out + last - first, work);
    return work;
}

using ordinary_double_matrix = ordinary_matrix<Rcpp::NumericVector>;

template <>
const int* ordinary_double_matrix::get_col(size_t c, int* work, size_t first, size_t last) {
    auto out = reader.get_col(c, first, last);
    std::copy(out, out + last - first, work);
    return work;
}

template <>
const double* ordinary_double_matrix::get_col(size_t c, double* work, size_t first, size_t last) {
    return reader.get_col(c, first, last);
}

/*** gCMatrix ***/

template <class V, typename TIT = typename V::iterator>
class gCMatrix : public sparse_lin_matrix {
public:
    gCMatrix(Rcpp::RObject mat) : reader(mat) {
        this->nrow = reader.get_nrow();
        this->ncol = reader.get_ncol();
        return;
    }
   
    ~gCMatrix() = default;
    gCMatrix(const gCMatrix&) = default;
    gCMatrix& operator=(const gCMatrix&) = default;
    gCMatrix(gCMatrix&&) = default;
    gCMatrix& operator=(gCMatrix&&) = default;

    const int* get_col(size_t c, int* work, size_t first, size_t last) {
        reader.get_col(c, work, first, last);
        return work;        
    }

    const int* get_row(size_t r, int* work, size_t first, size_t last) {
        reader.get_row(r, work, first, last);       
        return work;
    }

    const double* get_col(size_t c, double* work, size_t first, size_t last) {
        reader.get_col(c, work, first, last);
        return work;
    }

    const double* get_row(size_t r, double* work, size_t first, size_t last) {
        reader.get_row(r, work, first, last); 
        return work;
    }
    
    sparse_index<const int*, int> get_col(size_t c, int* work_x, int* work_i, size_t first, size_t last);

    sparse_index<const int*, int> get_row(size_t r, int* work_x, int* work_i, size_t first, size_t last) {
        return reader.template get_row<const int*>(r, work_x, work_i, first, last);
    }

    sparse_index<const double*, int> get_col(size_t c, double* work_x, int* work_i, size_t first, size_t last);

    sparse_index<const double*, int> get_row(size_t r, double* work_x, int* work_i, size_t first, size_t last) {
        return reader.template get_row<const double*>(r, work_x, work_i, first, last);
    }
private:
    gCMatrix_reader<V, TIT> reader;
};

using lgCMatrix = gCMatrix<Rcpp::LogicalVector, const int*>;

template <>
sparse_index<const int*, int> lgCMatrix::get_col(size_t c, int* work_x, int* work_i, size_t first, size_t last) {
    return reader.get_col(c, first, last);
}

template <>
sparse_index<const double*, int> lgCMatrix::get_col(size_t c, double* work_x, int* work_i, size_t first, size_t last) {
    return transplant<const double*>(reader.get_col(c, first, last), work_x, work_i);
}

using dgCMatrix = gCMatrix<Rcpp::NumericVector, const double*>;

template <>
sparse_index<const int*, int> dgCMatrix::get_col(size_t c, int* work_x, int* work_i, size_t first, size_t last) {
    return transplant<const int*>(reader.get_col(c, first, last), work_x, work_i);
}

template <>
sparse_index<const double*, int> dgCMatrix::get_col(size_t c, double* work_x, int* work_i, size_t first, size_t last) {
    return reader.get_col(c, first, last);
}

/*** SparseArraySeed ***/

template <class V, typename TIT = typename V::iterator>
class lin_SparseArraySeed : public sparse_lin_matrix {
public:
    lin_SparseArraySeed(Rcpp::RObject mat) : reader(mat) {
        this->nrow = reader.get_nrow();
        this->ncol = reader.get_ncol();
        return;
    }

    ~lin_SparseArraySeed() = default;
    lin_SparseArraySeed(const lin_SparseArraySeed&) = default;
    lin_SparseArraySeed& operator=(const lin_SparseArraySeed&) = default;
    lin_SparseArraySeed(lin_SparseArraySeed&&) = default;
    lin_SparseArraySeed& operator=(lin_SparseArraySeed&&) = default;

    const int* get_col(size_t c, int* work, size_t first, size_t last) {
        reader.get_col(c, work, first, last);
        return work;        
    }

    const int* get_row(size_t r, int* work, size_t first, size_t last) {
        reader.get_row(r, work, first, last);       
        return work;
    }

    const double* get_col(size_t c, double* work, size_t first, size_t last) {
        reader.get_col(c, work, first, last);
        return work;
    }

    const double* get_row(size_t r, double* work, size_t first, size_t last) {
        reader.get_row(r, work, first, last); 
        return work;
    }

    sparse_index<const int*, int> get_col(size_t c, int* work_x, int* work_i, size_t first, size_t last);

    sparse_index<const int*, int> get_row(size_t r, int* work_x, int* work_i, size_t first, size_t last) {
        return reader.template get_row<const int*>(r, work_x, work_i, first, last);
    }

    sparse_index<const double*, int> get_col(size_t c, double* work_x, int* work_i, size_t first, size_t last);

    sparse_index<const double*, int> get_row(size_t r, double* work_x, int* work_i, size_t first, size_t last) {
        return reader.template get_row<const double*>(r, work_x, work_i, first, last);
    }
private:
    SparseArraySeed_reader<V, TIT> reader;
};

using integer_SparseArraySeed = lin_SparseArraySeed<Rcpp::IntegerVector, const int*>;

template <>
sparse_index<const int*, int> integer_SparseArraySeed::get_col(size_t c, int* work_x, int* work_i, size_t first, size_t last) {
    return reader.get_col(c, first, last);
}

template <>
sparse_index<const double*, int> integer_SparseArraySeed::get_col(size_t c, double* work_x, int* work_i, size_t first, size_t last)
{
    return transplant<const double*>(reader.get_col(c, first, last), work_x, work_i);
}

using logical_SparseArraySeed = lin_SparseArraySeed<Rcpp::LogicalVector, const int*>;

template <>
sparse_index<const int*, int> logical_SparseArraySeed::get_col(size_t c, int* work_x, int* work_i, size_t first, size_t last) {
    return reader.get_col(c, first, last);
}

template <>
sparse_index<const double*, int> logical_SparseArraySeed::get_col(size_t c, double* work_x, int* work_i, size_t first, size_t last)
{
    return transplant<const double*>(reader.get_col(c, first, last), work_x, work_i);
}

using double_SparseArraySeed = lin_SparseArraySeed<Rcpp::NumericVector, const double*>;

template <>
sparse_index<const int*, int> double_SparseArraySeed::get_col(size_t c, int* work_x, int* work_i, size_t first, size_t last) {
    return transplant<const int*>(reader.get_col(c, first, last), work_x, work_i);
}

template <>
sparse_index<const double*, int> double_SparseArraySeed::get_col(size_t c, double* work_x, int* work_i, size_t first, size_t last)
{
    return reader.get_col(c, first, last);
}

}

#endif
