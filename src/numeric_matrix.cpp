#include "numeric_matrix.h"

namespace beachmat { 

/* Csparse numeric input methods. */

template<>
double Csparse_matrix<double, Rcpp::NumericVector>::get_empty() const { return 0; }

/* HDF5Matrix input methods. */

template<>
double HDF5_lin_matrix<double, Rcpp::NumericVector, REALSXP>::get(size_t r, size_t c) {
    double out;
    mat.extract_one(r, c, &out, H5::PredType::NATIVE_DOUBLE);
    return out; 
}

/* DelayedMatrix input methods. */

const std::vector<std::string> allowed_seeds={"dgeMatrix", "dgCMatrix", "dgTMatrix", "dspMatrix", "RleMatrix"};

template<>
std::unique_ptr<numeric_matrix> delayed_lin_matrix<double, Rcpp::NumericVector>::generate_seed(Rcpp::RObject incoming) {
    Rcpp::RObject seed=extract_seed(incoming, allowed_seeds);
    if (seed!=R_NilValue) {
        return create_numeric_matrix(seed);
    } else {
        return nullptr;
    }
} 

/* Sparse numeric output methods. */

template<>
double Csparse_output<double, Rcpp::NumericVector>::get_empty() const { return 0; }

/* HDF5 numeric output methods. */

template<>
double HDF5_output<double, Rcpp::NumericVector>::get_empty() const { return 0; }

template<>
Rcpp::RObject HDF5_output<double, Rcpp::NumericVector>::get_firstval() { 
    double first;
    extract_one(0, 0, &first);
    return Rcpp::NumericVector::create(first);
}

/* Dispatch definition */

std::unique_ptr<numeric_matrix> create_numeric_matrix(const Rcpp::RObject& incoming) { 
    if (incoming.isS4()) {
        std::string ctype=get_class(incoming);
        if (ctype=="dgeMatrix") { 
            return std::unique_ptr<numeric_matrix>(new dense_numeric_matrix(incoming));
        } else if (ctype=="dgCMatrix") { 
            return std::unique_ptr<numeric_matrix>(new Csparse_numeric_matrix(incoming));
        } else if (ctype=="dgTMatrix") {
            throw std::runtime_error("dgTMatrix not supported, convert to dgCMatrix");
        } else if (ctype=="dspMatrix") {
            return std::unique_ptr<numeric_matrix>(new Psymm_numeric_matrix(incoming));
        } else if (ctype=="HDF5Matrix") {
            return std::unique_ptr<numeric_matrix>(new HDF5_numeric_matrix(incoming));
        } else if (ctype=="RleMatrix") {
            return std::unique_ptr<numeric_matrix>(new Rle_numeric_matrix(incoming));
        } else if (ctype=="DelayedMatrix") { 
            return std::unique_ptr<numeric_matrix>(new delayed_numeric_matrix(incoming));
        }
        throw_custom_error("unsupported class '", ctype, "' for numeric_matrix");
    } 
    return std::unique_ptr<numeric_matrix>(new simple_numeric_matrix(incoming));
}

/* Output dispatch definition */

std::unique_ptr<numeric_output> create_numeric_output(int nrow, int ncol, const output_param& param) {
    switch (param.get_mode()) {
        case SIMPLE:
            return std::unique_ptr<numeric_output>(new simple_numeric_output(nrow, ncol));
        case SPARSE:
            return std::unique_ptr<numeric_output>(new sparse_numeric_output(nrow, ncol));
        case HDF5:
            return std::unique_ptr<numeric_output>(new HDF5_numeric_output(nrow, ncol, 
                        param.get_chunk_nrow(), param.get_chunk_ncol(), param.get_compression()));
        default:
            throw std::runtime_error("unsupported output mode for numeric matrices");
    }
}

}
