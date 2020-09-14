#ifndef BEACHMAT_BEACHMAT_H
#define BEACHMAT_BEACHMAT_H

/**
 * @file beachmat.h
 *
 * The main header file that should be `include`d to use the **beachmat** API.
 * This provides functions to easily construct an instance of the appropriate subclass from an R object.
 */

#include "Rcpp.h"
#include "lin_matrix.h"
#include <stdexcept>
#include <memory>

namespace beachmat {

/**
 * @internal
 *
 * Read a sparse logical, integer or numeric block into an instance of a `M` class.
 * `M` should most typically be either a `lin_matrix` or `sparse_lin_matrix`.
 *
 * @param block An R object containing a `dgCMatrix`, `lgCMatrix` or `SparseArraySeed`.
 *
 * @return A pointer to an instance of the `M` class.
 */
template <class M> 
std::unique_ptr<M> read_sparse_lin_block_raw (Rcpp::RObject block) {
    std::string ctype=get_class_name(block);

    if (ctype == "SparseArraySeed") {
        Rcpp::RObject nzdata=get_safe_slot(block, "nzdata");
        auto sexptype = nzdata.sexp_type();
        
        if (sexptype == INTSXP) {
            return std::unique_ptr<M>(new integer_SparseArraySeed(block));
        } else if (sexptype == REALSXP) {
            return std::unique_ptr<M>(new double_SparseArraySeed(block));
        } else if (sexptype == LGLSXP) {
            return std::unique_ptr<M>(new logical_SparseArraySeed(block));
        }

    } else if (ctype == "lgCMatrix") {
        return std::unique_ptr<M>(new lgCMatrix(block));

    } else if (ctype == "dgCMatrix") {
        return std::unique_ptr<M>(new dgCMatrix(block));

    } 

    return std::unique_ptr<M>();
}

/**
 * Read a logical, integer or numeric block into an instance of a `lin_matrix` subclass.
 * This can then be used to perform class- and type-agnostic extraction of row/column vectors.
 *
 * @param block An R object containing an ordinary submatrix, `dgCMatrix`, `lgCMatrix` or `SparseArraySeed`.
 *
 * @return A pointer to a `lin_matrix` instance.
 * This function will automatically choose the most appropriate subclass or throw an error if none are available.
 */
inline std::unique_ptr<lin_matrix> read_lin_block(Rcpp::RObject block) {
    if (block.isS4()) {
        auto ptr = read_sparse_lin_block_raw<lin_matrix>(block);
        if (ptr) {
            return ptr;
        }
    } else {
        auto sexptype = block.sexp_type();
        if (sexptype == INTSXP) {
            return std::unique_ptr<lin_matrix>(new ordinary_integer_matrix(block));
        } else if (sexptype == REALSXP) {
            return std::unique_ptr<lin_matrix>(new ordinary_double_matrix(block));
        } else if (sexptype == LGLSXP) {
            return std::unique_ptr<lin_matrix>(new ordinary_logical_matrix(block));
        }
    }

    throw std::runtime_error("'block' is not a recognized matrix representation");
}

/**
 * Read a sparse logical, integer or numeric block into an instance of a `sparse_lin_matrix` subclass.
 * This can then be used to perform class- and type-agnostic extraction of row/column vectors and their non-zero values.
 *
 * @param block An R object containing a `dgCMatrix`, `lgCMatrix` or `SparseArraySeed`.
 *
 * @return A pointer to a `sparse_lin_matrix` instance.
 * This function will automatically choose the most appropriate subclass or throw an error if none are available.
 */
inline std::unique_ptr<sparse_lin_matrix> read_sparse_lin_block(Rcpp::RObject block) {
    if (block.isS4()) {
        auto ptr = read_sparse_lin_block_raw<sparse_lin_matrix>(block);
        if (ptr) {
            return ptr;
        }
    }

    std::string ctype=get_class_name(block);
    throw std::runtime_error(ctype + std::string(" is not a recognized sparse representation"));
}

}

#endif
