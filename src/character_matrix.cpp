#include "character_matrix.h"

namespace beachmat {

/* Methods for the virtual class. */

character_matrix::character_matrix() {}

character_matrix::~character_matrix() {}

void character_matrix::get_col(size_t c, Rcpp::StringVector::iterator out) { 
    get_col(c, out, 0, get_nrow());
    return;
}

void character_matrix::get_row(size_t r, Rcpp::StringVector::iterator out) { 
    get_row(r, out, 0, get_ncol());
    return;
}

void character_matrix::get_cols(Rcpp::IntegerVector::iterator it, size_t n, Rcpp::StringVector::iterator out) {
    get_cols(it, n, out, 0, get_nrow());
    return;
}

void character_matrix::get_rows(Rcpp::IntegerVector::iterator it, size_t n, Rcpp::StringVector::iterator out) {
    get_rows(it, n, out, 0, get_ncol());
    return;
}

// Specialized getters

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
    return reader.get_const_col(c, first, last);
}

std::unique_ptr<character_matrix> simple_character_matrix::clone() const {
    return std::unique_ptr<character_matrix>(new simple_character_matrix(*this));
}

/* Methods for the helper class for the HDF5 character interface. */

HDF5_character_reader::HDF5_character_reader(const Rcpp::RObject& incoming) : HDF5_reader<Rcpp::String, STRSXP>(incoming), 
        str_type(this->get_datatype()) {

    if (str_type.isVariableStr()) { 
        throw std::runtime_error("variable-length strings not supported for HDF5_character_matrix");
    }

    bufsize=str_type.getSize(); 
    buffer.resize(bufsize * std::max({ this->get_ncol(), this->get_nrow(), size_t(1) }));
    return;
}

HDF5_character_reader::~HDF5_character_reader() {}

void HDF5_character_reader::get_row(size_t r, Rcpp::StringVector::iterator out, size_t first, size_t last) { 
    char* ref=buffer.data();
    this->extract_row(r, ref, str_type, first, last);
    for (size_t c=first; c<last; ++c, ref+=bufsize, ++out) {
        (*out)=ref; 
    }
    return;
} 

void HDF5_character_reader::get_col(size_t c, Rcpp::StringVector::iterator out, size_t first, size_t last) { 
    char* ref=buffer.data();
    this->extract_col(c, ref, str_type, first, last);
    for (size_t r=first; r<last; ++r, ref+=bufsize, ++out) {
        (*out)=ref; 
    }
    return;
}
 
Rcpp::String HDF5_character_reader::get(size_t r, size_t c) { 
    char* ref=buffer.data();
    this->extract_one(r, c, ref, str_type);
    return ref;
}

void HDF5_character_reader::get_rows(Rcpp::IntegerVector::iterator it, size_t n, Rcpp::StringVector::iterator out, size_t first, size_t last) {
    const size_t required=bufsize * n * (last - first);
    if (required > buffer.size()) {
        dim_checker::check_subset(first, last, get_ncol(), "column");
        buffer.resize(required);
    }

    char* ref=buffer.data();
    this->extract_rows(it, n, ref, str_type, first, last);
    for (size_t i=0; i<required; i+=bufsize, ref+=bufsize, ++out) {
        (*out)=ref; 
    }
    return;
}

void HDF5_character_reader::get_cols(Rcpp::IntegerVector::iterator it, size_t n, Rcpp::StringVector::iterator out, size_t first, size_t last) {
    const size_t required=bufsize * n * (last - first);
    if (required > buffer.size()) {
        dim_checker::check_subset(first, last, get_nrow(), "row");
        buffer.resize(required);
    }

    char* ref=buffer.data();
    this->extract_cols(it, n, ref, str_type, first, last);
    for (size_t i=0; i<required; i+=bufsize, ref+=bufsize, ++out) {
        (*out)=ref; 
    }
    return;
}

/* Methods for the Delayed character interface. */

std::unique_ptr<character_matrix> create_character_matrix_internal(const Rcpp::RObject&, bool); 

template<>
std::unique_ptr<character_matrix> delayed_character_reader::generate_seed(Rcpp::RObject incoming) {
    return create_character_matrix_internal(incoming, false);
}

/* Dispatch definition */

std::unique_ptr<character_matrix> create_character_matrix_internal(const Rcpp::RObject& incoming, bool delayed) { 
    if (incoming.isS4()) { 
        std::string ctype=get_class(incoming);
        if (ctype=="HDF5Matrix") {
            return std::unique_ptr<character_matrix>(new HDF5_character_matrix(incoming));
        } else if (delayed && ctype=="DelayedMatrix") { 
            return std::unique_ptr<character_matrix>(new delayed_character_matrix(incoming));
        }
        return std::unique_ptr<character_matrix>(new unknown_character_matrix(incoming));
    } 
    return std::unique_ptr<character_matrix>(new simple_character_matrix(incoming));
}

std::unique_ptr<character_matrix> create_character_matrix(const Rcpp::RObject& incoming) { 
    return create_character_matrix_internal(incoming, true);
}


}
