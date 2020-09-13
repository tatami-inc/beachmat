#ifndef BEACHMAT_BEACHMAT_H
#define BEACHMAT_BEACHMAT_H

#include "Rcpp.h"
#include "ordinary_matrix.h"
#include "Csparse_matrix.h"
#include <stdexcept>
#include <memory>

namespace beachmat {

/* read_block */

using integer_matrix = any_matrix<int>;

inline std::unique_ptr<integer_matrix> read_integer_block(Rcpp::RObject block) {
    if (block.isS4()) {
        std::string ctype=get_class_name(block);
        if (ctype=="SparseArraySeed") {
            return std::unique_ptr<integer_matrix>(new sparse_seed<Rcpp::IntegerVector>(block));
        }
        throw std::runtime_error(ctype + std::string(" is not a recognized integer matrix representation"));
    }
    return std::unique_ptr<integer_matrix>(new ordinary_matrix<Rcpp::IntegerVector>(block));
}

using double_matrix = any_matrix<double>;

inline std::unique_ptr<double_matrix> read_double_block(Rcpp::RObject block) {
    if (block.isS4()) {
        std::string ctype=get_class_name(block);
        if (ctype=="dgCMatrix") { 
            return std::unique_ptr<double_matrix>(new gCMatrix<Rcpp::NumericVector>(block));
        } else if (ctype=="SparseArraySeed") {
            return std::unique_ptr<double_matrix>(new sparse_seed<Rcpp::NumericVector>(block));
        }
        throw std::runtime_error(ctype + std::string(" is not a recognized numeric matrix representation"));
    }
    return std::unique_ptr<double_matrix>(new ordinary_matrix<Rcpp::NumericVector>(block));
}

using logical_matrix = any_matrix<int>;

inline std::unique_ptr<logical_matrix> read_logical_block(Rcpp::RObject block) {
    if (block.isS4()) {
        std::string ctype=get_class_name(block);
        if (ctype=="lgCMatrix") { 
            return std::unique_ptr<logical_matrix>(new gCMatrix<Rcpp::LogicalVector>(block));
        } else if (ctype=="SparseArraySeed") {
            return std::unique_ptr<logical_matrix>(new sparse_seed<Rcpp::LogicalVector>(block));
        }
        throw std::runtime_error(ctype + std::string(" is not a recognized logical matrix representation"));
    }
    return std::unique_ptr<logical_matrix>(new ordinary_matrix<Rcpp::LogicalVector>(block));
}

/* read_sparse_block */

using integer_Csparse_matrix = Csparse_matrix<int>;

inline std::unique_ptr<integer_Csparse_matrix> read_integer_Csparse_block(Rcpp::RObject block) {
    std::string ctype=get_class_name(block);
    if (ctype=="SparseArraySeed") {
        return std::unique_ptr<integer_Csparse_matrix>(new sparse_seed<Rcpp::IntegerVector>(block));
    }
    throw std::runtime_error(ctype + std::string(" is not a recognized integer matrix representation"));
}

using double_Csparse_matrix = Csparse_matrix<double>;

inline std::unique_ptr<double_Csparse_matrix> read_double_Csparse_block(Rcpp::RObject block) {
    std::string ctype=get_class_name(block);
    if (ctype=="dgCMatrix") { 
        return std::unique_ptr<double_Csparse_matrix>(new gCMatrix<Rcpp::NumericVector>(block));
    } else if (ctype=="SparseArraySeed") {
        return std::unique_ptr<double_Csparse_matrix>(new sparse_seed<Rcpp::NumericVector>(block));
    }
    throw std::runtime_error(ctype + std::string(" is not a recognized numeric matrix representation"));
}

using logical_Csparse_matrix = Csparse_matrix<int>;

inline std::unique_ptr<logical_Csparse_matrix> read_logical_Csparse_block(Rcpp::RObject block) {
    std::string ctype=get_class_name(block);
    if (ctype=="lgCMatrix") { 
        return std::unique_ptr<logical_Csparse_matrix>(new gCMatrix<Rcpp::LogicalVector>(block));
    } else if (ctype=="SparseArraySeed") {
        return std::unique_ptr<logical_Csparse_matrix>(new sparse_seed<Rcpp::LogicalVector>(block));
    }
    throw std::runtime_error(ctype + std::string(" is not a recognized logical matrix representation"));
}

}

#endif
