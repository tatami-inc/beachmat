#include "character_matrix.h"

namespace beachmat {

/* Methods for the virtual class. */

character_matrix::character_matrix() {}

character_matrix::~character_matrix() {}

void character_matrix::get_col(size_t c, Rcpp::StringVector::iterator out) { 
    get_col(c, out, 0, get_nrow());
}

void character_matrix::get_row(size_t r, Rcpp::StringVector::iterator out) { 
    get_row(r, out, 0, get_ncol());
}

Rcpp::StringVector::iterator character_matrix::get_const_col(size_t c, Rcpp::StringVector::iterator work) {
    return get_const_col(c, work, 0, get_nrow());
}

Rcpp::StringVector::iterator character_matrix::get_const_col(size_t c, Rcpp::StringVector::iterator work, size_t first, size_t last) {
    get_col(c, work, first, last);
    return work;
}

const_col_indexed_info<Rcpp::StringVector> character_matrix::get_const_col_indexed(size_t c, Rcpp::StringVector::iterator work) {
    return get_const_col_indexed(c, work, 0, get_nrow());
}

const_col_indexed_info<Rcpp::StringVector> character_matrix::get_const_col_indexed(size_t c, Rcpp::StringVector::iterator work, size_t first, size_t last) {
    if (indices.size()!=this->get_nrow()) {
        indices=Rcpp::IntegerVector(this->get_nrow());
        std::iota(indices.begin(), indices.end(), 0); // populating with indices.
    }
    return const_col_indexed_info<Rcpp::StringVector>(last - first, indices.begin() + first, get_const_col(c, work, first, last));
}

/* Methods for the simple character matrix. */

simple_character_matrix::simple_character_matrix(const Rcpp::RObject& incoming) : simple_character_precursor (incoming) {}

simple_character_matrix::~simple_character_matrix() {}

Rcpp::StringVector::iterator simple_character_matrix::get_const_col(size_t c, Rcpp::StringVector::iterator work, size_t first, size_t last) {
    return mat.get_const_col(c, first, last);
}

std::unique_ptr<character_matrix> simple_character_matrix::clone() const {
    return std::unique_ptr<character_matrix>(new simple_character_matrix(*this));
}

/* Methods for the HDF5 character matrix. */

HDF5_character_matrix::HDF5_character_matrix(const Rcpp::RObject& incoming) : mat(incoming), str_type(mat.get_datatype()) {
    if (str_type.isVariableStr()) { 
        throw std::runtime_error("variable-length strings not supported for HDF5_character_matrix");
    }
    bufsize=str_type.getSize(); 
    row_buf.resize(bufsize*(mat.get_ncol()));
    col_buf.resize(bufsize*(mat.get_nrow()));
    one_buf.resize(bufsize);
    return;
}

HDF5_character_matrix::~HDF5_character_matrix() {}

size_t HDF5_character_matrix::get_nrow() const {
    return mat.get_nrow();
}

size_t HDF5_character_matrix::get_ncol() const {
    return mat.get_ncol();
}

void HDF5_character_matrix::get_row(size_t r, Rcpp::StringVector::iterator out, size_t first, size_t last) { 
    char* ref=row_buf.data();
    mat.extract_row(r, ref, str_type, first, last);
    for (size_t c=first; c<last; ++c, ref+=bufsize, ++out) {
        (*out)=ref; 
    }
    return;
} 

void HDF5_character_matrix::get_col(size_t c, Rcpp::StringVector::iterator out, size_t first, size_t last) { 
    char* ref=col_buf.data();
    mat.extract_col(c, ref, str_type, first, last);
    for (size_t r=first; r<last; ++r, ref+=bufsize, ++out) {
        (*out)=ref; 
    }
    return;
}
 
Rcpp::String HDF5_character_matrix::get(size_t r, size_t c) { 
    char* ref=one_buf.data();
    mat.extract_one(r, c, ref, str_type);
    return ref;
}

std::unique_ptr<character_matrix> HDF5_character_matrix::clone() const {
    return std::unique_ptr<character_matrix>(new HDF5_character_matrix(*this));
}
 
Rcpp::RObject HDF5_character_matrix::yield() const {
    return mat.yield();
}

matrix_type HDF5_character_matrix::get_matrix_type() const {
    return mat.get_matrix_type();
}

/* Methods for the delayed character matrix. */

const std::vector<std::string> allowed_s4_seeds={"RleMatrix"};

template<>
std::unique_ptr<character_matrix> delayed_character_helper::generate_seed(Rcpp::RObject incoming) {
    Rcpp::RObject seed=extract_seed(incoming, allowed_s4_seeds);
    if (seed!=R_NilValue) {
        return create_character_matrix(seed);
    } else {
        return nullptr;
    }
}

template <>
std::unique_ptr<character_matrix> delayed_character_helper::generate_unknown_seed(Rcpp::RObject incoming) {
    return std::unique_ptr<character_matrix>(new unknown_character_matrix(incoming));
}

/* Dispatch definition */

std::unique_ptr<character_matrix> create_character_matrix(const Rcpp::RObject& incoming) { 
    if (incoming.isS4()) { 
        std::string ctype=get_class(incoming);
        if (ctype=="HDF5Matrix") {
            return std::unique_ptr<character_matrix>(new HDF5_character_matrix(incoming));
        } else if (ctype=="RleMatrix") { 
            return std::unique_ptr<character_matrix>(new Rle_character_matrix(incoming));
        } else if (ctype=="DelayedMatrix") { 
            return std::unique_ptr<character_matrix>(new delayed_character_matrix(incoming));
        }
        return std::unique_ptr<character_matrix>(new unknown_character_matrix(incoming));
    } 
    return std::unique_ptr<character_matrix>(new simple_character_matrix(incoming));
}

}
