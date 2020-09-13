#ifndef BEACHMAT_BEACHMAT_H
#define BEACHMAT_BEACHMAT_H

#include "Rcpp.h"
#include "ordinary_matrix.h"
#include "Csparse_matrix.h"
#include <stdexcept>

namespace beachmat {

/* read_block */

typename any_matrix<Rcpp::IntegerVector> integer_matrix;

std::unique_ptr<integer_matrix*> read_integer_block(Rcpp::RObject block) {
    if (block.isS4()) {
        std::string ctype=get_class_name(block);
        if (ctype=="SparseArraySeed") {
            return std::unique_ptr<integer_matrix*>(new sparse_seed<int, Rcpp::IntegerVector>(block));
        }
        throw std::runtime_error(ctype + std::string(" is not a recognized integer matrix representation"));
    }
    return std::unique_ptr<integer_matrix*>(new ordinary_matrix<int, Rcpp::IntegerVector>(block));
}

typename any_matrix<Rcpp::NumericVector> double_matrix;

std::unique_ptr<numeric_matrix*> read_double_block(Rcpp::RObject block) {
    if (block.isS4()) {
        std::string ctype=get_class_name(block);
        if (ctype=="dgCMatrix") { 
            return std::unique_ptr<double_matrix*>(new gCMatrix<double, Rcpp::NumericVector>(block));
        } else if (ctype=="SparseArraySeed") {
            return std::unique_ptr<double_matrix*>(new sparse_seed<double, Rcpp::NumericVector>(block));
        }
        throw std::runtime_error(ctype + std::string(" is not a recognized numeric matrix representation"));
    }
    return std::unique_ptr<double_matrix*>(new ordinary_matrix<double, Rcpp::NumericVector>(block));
}

typename any_matrix<Rcpp::LogicalVector> logical_matrix;

std::unique_ptr<logical_matrix*> read_logical_block(Rcpp::RObject block) {
    if (block.isS4()) {
        std::string ctype=get_class_name(block);
        if (ctype=="lgCMatrix") { 
            return std::unique_ptr<logical_matrix*>(new gCMatrix<int, Rcpp::LogicalVector>(block));
        } else if (ctype=="SparseArraySeed") {
            return std::unique_ptr<logical_matrix*>(new sparse_seed<int, Rcpp::LogicalVector>(block));
        }
        throw std::runtime_error(ctype + std::string(" is not a recognized logical matrix representation"));
    }
    return std::unique_ptr<logical_matrix*>(new ordinary_matrix<logical, Rcpp::LogicalVector>(block));
}

/* read_sparse_block */

typename Csparse_matrix<Rcpp::IntegerVector> integer_Csparse_matrix;

std::unique_ptr<integer_Csparse_matrix*> read_integer_Csparse_block(Rcpp::RObject block) {
    std::string ctype=get_class_name(block);
    if (ctype=="SparseArraySeed") {
        return std::unique_ptr<integer_Csparse_matrix*>(new sparse_seed<int, Rcpp::IntegerVector>(block));
    }
    throw std::runtime_error(ctype + std::string(" is not a recognized integer matrix representation"));
}

typename Csparse_matrix<Rcpp::NumericVector> double_Csparse_matrix;

std::unique_ptr<numeric_Csparse_matrix*> read_double_Csparse_block(Rcpp::RObject block) {
    std::string ctype=get_class_name(block);
    if (ctype=="dgCMatrix") { 
        return std::unique_ptr<double_Csparse_matrix*>(new gCMatrix<double, Rcpp::NumericVector>(block));
    } else if (ctype=="SparseArraySeed") {
        return std::unique_ptr<double_Csparse_matrix*>(new sparse_seed<double, Rcpp::NumericVector>(block));
    }
    throw std::runtime_error(ctype + std::string(" is not a recognized numeric matrix representation"));
}

typename Csparse_matrix<Rcpp::LogicalVector> logical_Csparse_matrix;

std::unique_ptr<logical_Csparse_matrix*> read_logical_Csparse_block(Rcpp::RObject block) {
    std::string ctype=get_class_name(block);
    if (ctype=="lgCMatrix") { 
        return std::unique_ptr<logical_Csparse_matrix*>(new gCMatrix<int, Rcpp::LogicalVector>(block));
    } else if (ctype=="SparseArraySeed") {
        return std::unique_ptr<logical_Csparse_matrix*>(new sparse_seed<int, Rcpp::LogicalVector>(block));
    }
    throw std::runtime_error(ctype + std::string(" is not a recognized logical matrix representation"));
}

}
