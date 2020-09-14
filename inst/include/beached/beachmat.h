#ifndef BEACHMAT_BEACHMAT_H
#define BEACHMAT_BEACHMAT_H

#include "Rcpp.h"
#include "lin_matrix.h"
#include <stdexcept>
#include <memory>

namespace beachmat {

/* read_block */

inline std::unique_ptr<lin_matrix> read_lin_block(Rcpp::RObject block) {
    if (block.isS4()) {
        std::string ctype=get_class_name(block);

        if (ctype == "SparseArraySeed") {
            Rcpp::RObject nzdata=get_safe_slot(block, "nzdata");
            auto sexptype = nzdata.sexp_type();
            
            if (sexptype == INTSXP) {
                return std::unique_ptr<lin_matrix>(new integer_SparseArraySeed(block));
            } else if (sexptype == REALSXP) {
                return std::unique_ptr<lin_matrix>(new double_SparseArraySeed(block));
            } else if (sexptype == LGLSXP) {
                return std::unique_ptr<lin_matrix>(new logical_SparseArraySeed(block));
            }

        } else if (ctype == "lgCMatrix") {
            return std::unique_ptr<lin_matrix>(new lgCMatrix(block));

        } else if (ctype == "dgCMatrix") {
            return std::unique_ptr<lin_matrix>(new dgCMatrix(block));

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

inline std::unique_ptr<sparse_lin_matrix> read_sparse_lin_block(Rcpp::RObject block) {
    std::string ctype=get_class_name(block);

    if (ctype == "SparseArraySeed") {
        Rcpp::RObject nzdata=get_safe_slot(block, "nzdata");
        auto sexptype = nzdata.sexp_type();

        if (sexptype == INTSXP) {
            return std::unique_ptr<sparse_lin_matrix>(new integer_SparseArraySeed(block));
        } else if (sexptype == REALSXP) {
            return std::unique_ptr<sparse_lin_matrix>(new double_SparseArraySeed(block));
        } else if (sexptype == LGLSXP) {
            return std::unique_ptr<sparse_lin_matrix>(new logical_SparseArraySeed(block));
        }

    } else if (ctype == "lgCMatrix") {
        return std::unique_ptr<sparse_lin_matrix>(new lgCMatrix(block));

    } else if (ctype == "dgCMatrix") {
        return std::unique_ptr<sparse_lin_matrix>(new dgCMatrix(block));

    }

    throw std::runtime_error(ctype + std::string(" is not a recognized sparse representation"));
}

}

#endif
