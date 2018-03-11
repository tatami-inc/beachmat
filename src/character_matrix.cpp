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

character_matrix::const_col_indexed_info 
character_matrix::get_const_col_indexed(size_t c, Rcpp::StringVector::iterator work) {
    return get_const_col_indexed(c, work, 0, get_nrow());
}

character_matrix::const_col_indexed_info 
character_matrix::get_const_col_indexed(size_t c, Rcpp::StringVector::iterator work, size_t first, size_t last) {
    if (indices.size()!=this->get_nrow()) {
        indices=Rcpp::IntegerVector(this->get_nrow());
        std::iota(indices.begin(), indices.end(), 0); // populating with indices.
    }
    return const_col_indexed_info(last - first, indices.begin() + first, get_const_col(c, work, first, last));
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

delayed_character_matrix::delayed_character_matrix(const Rcpp::RObject& incoming) : original(incoming), seed_ptr(nullptr) {
    check_DelayedMatrix(incoming);
    
    // Trying to generate the seed, if it's a valid object in itself.
    if (only_delayed_coord_changes(incoming)) {
        seed_ptr=generate_seed(incoming);
    }
        
    // If the seed is still NULL, we switch to a chunked matrix format.
    if (seed_ptr.get()==NULL) { 
        seed_ptr=std::unique_ptr<character_matrix>(new delayed_character_matrix::enslaved(incoming));
    } else {
        transformer=delayed_coord_transformer<Rcpp::String, Rcpp::StringVector>(incoming, seed_ptr.get());
    }

    return;
}

delayed_character_matrix::~delayed_character_matrix() {}

delayed_character_matrix::delayed_character_matrix(const delayed_character_matrix& other) : original(other.original),
    seed_ptr(other.seed_ptr->clone()), transformer(other.transformer) {}

delayed_character_matrix& delayed_character_matrix::operator=(const delayed_character_matrix& other) {
    original=other.original;
    seed_ptr=other.seed_ptr->clone();
    transformer=other.transformer;
    return *this;
}

size_t delayed_character_matrix::get_nrow() const {
    return transformer.get_nrow(); 
}

size_t delayed_character_matrix::get_ncol() const {
    return transformer.get_ncol();
}

void delayed_character_matrix::get_col(size_t c, Rcpp::StringVector::iterator out, size_t first, size_t last) {
    transformer.get_col(seed_ptr.get(), c, out, first, last);
    return;
}

void delayed_character_matrix::get_row(size_t r, Rcpp::StringVector::iterator out, size_t first, size_t last) {
    transformer.get_row(seed_ptr.get(), r, out, first, last);
    return;
}

Rcpp::String delayed_character_matrix::get(size_t r, size_t c) {
    return transformer.get(seed_ptr.get(), r, c);
}

std::unique_ptr<character_matrix> delayed_character_matrix::clone() const {
    return std::unique_ptr<character_matrix>(new delayed_character_matrix(*this));
}

Rcpp::RObject delayed_character_matrix::yield() const {
    return original;
}

matrix_type delayed_character_matrix::get_matrix_type() const { // returns the type of the SEED!
    return DELAYED;
}

const std::vector<std::string> allowed_seeds={"RleMatrix"};

std::unique_ptr<character_matrix> delayed_character_matrix::generate_seed(Rcpp::RObject incoming) {
    Rcpp::RObject seed=extract_seed(incoming, allowed_seeds);
    if (seed!=R_NilValue) {
        return create_character_matrix(seed);
    } else {
        return nullptr;
    }
}

delayed_character_matrix::enslaved::enslaved(const Rcpp::RObject& in) : delayed_character_matrix::enslaved_precursor(in) {}

delayed_character_matrix::enslaved::~enslaved() {}

//template<typename T, class V>
//typename V::iterator delayed_character_matrix::enslaved::get_const_col(size_t c, typename V::iterator out, size_t first, size_t last) {
//    return (this->mat).get_const_col(c, out, first, last);
//}

std::unique_ptr<character_matrix > delayed_character_matrix::enslaved::clone() const {
    return std::unique_ptr<character_matrix>(new delayed_character_matrix::enslaved(*this));
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
        std::stringstream err;
        err << "unsupported class '" << ctype << "' for character_matrix";
        throw std::runtime_error(err.str().c_str());
    } 
    return std::unique_ptr<character_matrix>(new simple_character_matrix(incoming));
}

}
