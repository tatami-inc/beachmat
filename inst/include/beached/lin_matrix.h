#ifndef BEACHMAT_LIN_MATRIX_H
#define BEACHMAT_LIN_MATRIX_H

/**
 * @file lin_matrix.h
 *
 * Class definitions for the logical, integer or numeric (LIN) matrix.
 */

#include "Rcpp.h"
#include "ordinary_reader.h"
#include "Csparse_reader.h"
#include "utils.h"

namespace beachmat {

/**
 * @brief Virtual base class for a logical, integer or numeric (i.e., double-precision) matrix,
 * providing methods to extract rows or columns in dense form.
 */
class lin_matrix {
public:
    /**
     * Trivial constructor.
     */
    lin_matrix() {}

    virtual ~lin_matrix() = default;
    lin_matrix(const lin_matrix&) = default;
    lin_matrix& operator=(const lin_matrix&) = default;
    lin_matrix(lin_matrix&&) = default;
    lin_matrix& operator=(lin_matrix&&) = default;

    /**
     * Extract a column as an array of integers.
     *
     * @param c Index of the column of interest.
     * @param work The workspace, potentially used to store extracted values.
     * This should have at least `last - first` addressable elements.
     * @param first The index of the first row of interest.
     * @param last The index of the first row of that is _not_ of interest.
     *
     * @return A pointer is returned to the values of `c` as integers, starting at the `first` element.
     * This may or may not involve populating `work` with values copied from the underlying matrix;
     * a copy will have been performed if the return value compares equal to `work`.
     */
    virtual const int* get_col(size_t c, int* work, size_t first, size_t last) = 0;

    /**
     * Extract a row as an array of integers.
     *
     * @param r Index of the row of interest.
     * @param work The workspace, potentially used to store extracted values.
     * This should have at least `last - first` addressable elements.
     * @param first The index of the first column of interest.
     * @param last The index of the first column of that is _not_ of interest.
     *
     * @return A pointer is returned to the values of `r` as integers, starting at the `first` element.
     * This involves creating a copy in the workspace so the return value is always equal to `work`.
     */
    virtual const int* get_row(size_t r, int* work, size_t first, size_t last) = 0;

    /**
     * Extract a column as an array of doubles.
     *
     * @param c Index of the column of interest.
     * @param work The workspace, potentially used to store extracted values.
     * This should have at least `last - first` addressable elements.
     * @param first The index of the first row of interest.
     * @param last The index of the first row of that is _not_ of interest.
     *
     * @return A pointer is returned to the values of `c` as doubles, starting at the `first` element.
     * This may or may not involve populating `work` with values copied from the underlying matrix;
     * a copy will have been performed if the return value compares equal to `work`.
     */
    virtual const double* get_col(size_t c, double* work, size_t first, size_t last) = 0;

    /**
     * Extract a row as an array of doubles.
     *
     * @param r Index of the row of interest.
     * @param work The workspace, potentially used to store extracted values.
     * This should have at least `last - first` addressable elements.
     * @param first The index of the first column of interest.
     * @param last The index of the first column of that is _not_ of interest.
     *
     * @return A pointer is returned to the values of `r` as doubles, starting at the `first` element.
     * This involves creating a copy in the workspace so the return value is always equal to `work`.
     */
    virtual const double* get_row(size_t r, double* work, size_t first, size_t last) = 0;

    /**
     * Extract a column as an array of integers.
     *
     * @param c Index of the column of interest.
     * @param work The workspace, potentially used to store extracted values.
     * This should have at least `last - first` addressable elements.
     *
     * @return A pointer is returned to the values of `c` as integers, starting at the first element.
     * This may or may not involve populating `work` with values copied from the underlying matrix;
     * a copy will have been performed if the return value compares equal to `work`.
     */
    const int* get_col(size_t r, int* work) {
        return get_col(r, work, 0, nrow);
    }

    /**
     * Extract a row as an array of integers.
     *
     * @param r Index of the column of interest.
     * @param work The workspace, potentially used to store extracted values.
     * This should have at least `last - first` addressable elements.
     *
     * @return A pointer is returned to the values of `r` as integers, starting at the first element.
     * This involves creating a copy in the workspace so the return value is always equal to `work`.
     */
    const int* get_row(size_t r, int* work) {
        return get_row(r, work, 0, ncol);        
    }

    /**
     * Extract a column as an array of doubles.
     *
     * @param c Index of the column of interest.
     * @param work The workspace, potentially used to store extracted values.
     * This should have at least `last - first` addressable elements.
     *
     * @return A pointer is returned to the values of `c` as doubles, starting at the first element.
     * This may or may not involve populating `work` with values copied from the underlying matrix;
     * a copy will have been performed if the return value compares equal to `work`.
     */
    const double* get_col(size_t c, double* work) {
        return get_col(c, work, 0, nrow);
    }

    /**
     * Extract a row as an array of doubles.
     *
     * @param r Index of the column of interest.
     * @param work The workspace, potentially used to store extracted values.
     * This should have at least `last - first` addressable elements.
     *
     * @return A pointer is returned to the values of `r` as doubles, starting at the first element.
     * This involves creating a copy in the workspace so the return value is always equal to `work`.
     */
    const double* get_row(size_t r, double* work) {
        return get_row(r, work, 0, ncol);        
    }

    /**
     * Get the number of rows in the matrix.
     */
    size_t get_nrow() const { return nrow; }

    /**
     * Get the number of columns in the matrix.
     */
    size_t get_ncol() const { return ncol; }

    /**
     * Is the matrix sparse?
     */
    virtual bool is_sparse() const { return false; }
protected:
    size_t nrow=0, ncol=0;
};

/**
 * @brief Virtual base class for a sparse logical, integer or numeric (i.e., double-precision) matrix,
 * providing methods to extract rows or columns in dense or sparse form.
 */
class sparse_lin_matrix : public lin_matrix {
public:
    /**
     * Trivial constructor.
     */
    sparse_lin_matrix() {}

    ~sparse_lin_matrix() = default;
    sparse_lin_matrix(const sparse_lin_matrix&) = default;
    sparse_lin_matrix& operator=(const sparse_lin_matrix&) = default;
    sparse_lin_matrix(sparse_lin_matrix&&) = default;
    sparse_lin_matrix& operator=(sparse_lin_matrix&&) = default;

    /**
     * Extract all non-zero elements in a row, storing values as integers.
     *
     * @param r Index of the row of interest.
     * @param work_x The workspace for extracted non-zero values.
     * This should have at least `last - first` addressable elements.
     * @param work_i The workspace for column indices.
     * This should have at least `last - first` addressable elements.
     * @param first The index of the first column of interest.
     * @param last The index of the first column of that is _not_ of interest.
     *
     * @return A `sparse_index` is returned containing pointers to the workspace,
     * containing non-zero elements in `r` with column indices in `[first, last)`.
     * A copy of non-zero values and their indices is always performed into the workspace.
     */
    virtual sparse_index<const int*, int> get_row(size_t, int*, int*, size_t, size_t) = 0;

    /**
     * Extract all non-zero elements in a column, storing values as integers.
     *
     * @param r Index of the column of interest.
     * @param work_x The workspace for extracted non-zero values.
     * This should have at least `last - first` addressable elements.
     * @param work_i The workspace for column indices.
     * This should have at least `last - first` addressable elements.
     * @param first The index of the first row of interest.
     * @param last The index of the first row of that is _not_ of interest.
     *
     * @return A `sparse_index` is returned containing pointers to non-zero elements in `c` with row indices in `[first, last)`.
     * This may or may not involve populating `work` with values copied from the underlying matrix;
     * a copy will have been performed if the return value's pointers compare equal to `work_x` and `work_i`.
     */
    virtual sparse_index<const int*, int> get_col(size_t, int*, int*, size_t, size_t) = 0;

    /**
     * Extract all non-zero elements in a row, storing values as doubles.
     *
     * @param r Index of the row of interest.
     * @param work_x The workspace for extracted non-zero values.
     * This should have at least `last - first` addressable elements.
     * @param work_i The workspace for column indices.
     * This should have at least `last - first` addressable elements.
     * @param first The index of the first column of interest.
     * @param last The index of the first column of that is _not_ of interest.
     *
     * @return A `sparse_index` is returned containing pointers to the workspace.
     * containing non-zero elements in `r` with column indices in `[first, last)`.
     * A copy of non-zero values and their indices is always performed into the workspace.
     */
    virtual sparse_index<const double*, int> get_row(size_t, double*, int*, size_t, size_t) = 0;

    /**
     * Extract all non-zero elements in a column, storing values as doubles.
     *
     * @param r Index of the column of interest.
     * @param work_x The workspace for extracted non-zero values.
     * This should have at least `last - first` addressable elements.
     * @param work_i The workspace for column indices.
     * This should have at least `last - first` addressable elements.
     * @param first The index of the first row of interest.
     * @param last The index of the first row of that is _not_ of interest.
     *
     * @return A `sparse_index` is returned containing pointers to non-zero elements in `c` with row indices in `[first, last)`.
     * This may or may not involve populating `work` with values copied from the underlying matrix;
     * a copy will have been performed if the return value's pointers compare equal to `work_x` and `work_i`.
     */
    virtual sparse_index<const double*, int> get_col(size_t, double*, int*, size_t, size_t) = 0;

    /**
     * Extract all non-zero elements in a column, storing values as integers.
     *
     * @param r Index of the column of interest.
     * @param work_x The workspace for extracted non-zero values.
     * This should have at least `last - first` addressable elements.
     * @param work_i The workspace for column indices.
     * This should have at least `last - first` addressable elements.
     *
     * @return A `sparse_index` is returned containing pointers to all non-zero elements in `c`.
     * This may or may not involve populating `work` with values copied from the underlying matrix;
     * a copy will have been performed if the return value's pointers compare equal to `work_x` and `work_i`.
     */
    sparse_index<const int*, int> get_col(size_t c, int* work_x, int* work_i) {
        return get_col(c, work_x, work_i, 0, this->nrow);
    }

    /**
     * Extract all non-zero elements in a row, storing values as integers.
     *
     * @param r Index of the row of interest.
     * @param work_x The workspace for extracted non-zero values.
     * This should have at least `last - first` addressable elements.
     * @param work_i The workspace for column indices.
     * This should have at least `last - first` addressable elements.
     *
     * @return A `sparse_index` is returned containing pointers to the workspace,
     * containing all non-zero elements in `r`.
     * A copy of non-zero values and their indices is always performed into the workspace.
     */
    sparse_index<const int*, int> get_row(size_t r, int* work_x, int* work_i) {
        return get_row(r, work_x, work_i, 0, this->ncol);
    }

    /**
     * Extract all non-zero elements in a column, storing values as doubles.
     *
     * @param r Index of the column of interest.
     * @param work_x The workspace for extracted non-zero values.
     * This should have at least `last - first` addressable elements.
     * @param work_i The workspace for column indices.
     * This should have at least `last - first` addressable elements.
     *
     * @return A `sparse_index` is returned containing pointers to all non-zero elements in `c`.
     * This may or may not involve populating `work` with values copied from the underlying matrix;
     * a copy will have been performed if the return value's pointers compare equal to `work_x` and `work_i`.
     */
    sparse_index<const double*, int> get_col(size_t c, double* work_x, int* work_i) {
        return get_col(c, work_x, work_i, 0, this->nrow);
    }

    /**
     * Extract all non-zero elements in a row, storing values as doubles.
     *
     * @param r Index of the row of interest.
     * @param work_x The workspace for extracted non-zero values.
     * This should have at least `last - first` addressable elements.
     * @param work_i The workspace for column indices.
     * This should have at least `last - first` addressable elements.
     *
     * @return A `sparse_index` is returned containing pointers to the workspace,
     * containing all non-zero elements in `r`.
     * A copy of non-zero values and their indices is always performed into the workspace.
     */
    sparse_index<const double*, int> get_row(size_t r, double* work_x, int* work_i) {
        return get_row(r, work_x, work_i, 0, this->ncol);
    }

    bool is_sparse() const { return true; }
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
