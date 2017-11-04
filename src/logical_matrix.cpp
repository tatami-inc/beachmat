#include "logical_matrix.h"

namespace beachmat {

/* Csparse logical input methods. */

template<>
int Csparse_matrix<int, Rcpp::LogicalVector>::get_empty() const { return 0; }

/* DelayedMatrix input methods. */

template<>
std::unique_ptr<logical_matrix> delayed_lin_matrix<int, Rcpp::LogicalVector>::generate_seed(Rcpp::RObject incoming) {
    // incoming should be a "DelayedMatrix" object, not the seed within it!
    bool isokay=false;
    Rcpp::RObject seed(get_safe_slot(incoming, "seed"));

    if (seed.isS4()) { 
        std::string ctype=get_class(seed);
        if (ctype=="lgeMatrix"
                || ctype=="lgCMatrix" 
                || ctype=="lgTMatrix" 
                || ctype=="lspMatrix" 
                || ctype=="RleMatrix") {
            isokay=true;
            incoming=seed;
        } else if (ctype=="HDF5ArraySeed") {
            isokay=true;
            incoming=delayed_seed_to_HDF5Matrix(seed);
        }
    } else {
        isokay=true;
        incoming=seed;
    }

    if (isokay) {
        return create_logical_matrix(incoming);
    } else {
        return nullptr;
    }
} 

/* Sparse logical output methods. */

template<>
int Csparse_output<int, Rcpp::LogicalVector>::get_empty() const { return 0; }

/* HDF5 logical output methods. */

template<>
int HDF5_output<int, LGLSXP>::get_empty() const { return 0; }

template<>
Rcpp::RObject HDF5_output<int, LGLSXP>::get_firstval() { 
    int first;
    extract_one(0, 0, &first);
    return Rcpp::LogicalVector::create(first);
}

/* Dispatch definition */

std::unique_ptr<logical_matrix> create_logical_matrix(const Rcpp::RObject& incoming) { 
    if (incoming.isS4()) {
        std::string ctype=get_class(incoming);
        if (ctype=="lgeMatrix") { 
            return std::unique_ptr<logical_matrix>(new dense_logical_matrix(incoming));
        } else if (ctype=="lgCMatrix") { 
            return std::unique_ptr<logical_matrix>(new Csparse_logical_matrix(incoming));
        } else if (ctype=="lgTMatrix") {
            throw std::runtime_error("lgTMatrix not supported, convert to lgCMatrix");
        } else if (ctype=="lspMatrix") {
            return std::unique_ptr<logical_matrix>(new Psymm_logical_matrix(incoming));
        } else if (ctype=="HDF5Matrix") {
            return std::unique_ptr<logical_matrix>(new HDF5_logical_matrix(incoming));
        } else if (ctype=="RleMatrix") {
            return std::unique_ptr<logical_matrix>(new Rle_logical_matrix(incoming));
        } else if (ctype=="DelayedMatrix") { 
            return std::unique_ptr<logical_matrix>(new delayed_logical_matrix(incoming));
        }
        throw_custom_error("unsupported class '", ctype, "' for logical_matrix");
    } 
    return std::unique_ptr<logical_matrix>(new simple_logical_matrix(incoming));
}

/* Output dispatch definition */

std::unique_ptr<logical_output> create_logical_output(int nrow, int ncol, const output_param& param) {
    switch (param.get_mode()) {
        case SIMPLE:
            return std::unique_ptr<logical_output>(new simple_logical_output(nrow, ncol));
        case SPARSE:
            return std::unique_ptr<logical_output>(new sparse_logical_output(nrow, ncol));
        case HDF5:
            return std::unique_ptr<logical_output>(new HDF5_logical_output(nrow, ncol,
                        param.get_chunk_nrow(), param.get_chunk_ncol(), param.get_compression()));
        default:
            throw std::runtime_error("unsupported output mode for logical matrices");
    }
}

}
