#include "character_output.h"

namespace beachmat {

/* Methods for the virtual class. */

character_output::character_output() {}

character_output::~character_output() {}

void character_output::get_col(size_t c, Rcpp::StringVector::iterator out) { 
    get_col(c, out, 0, get_nrow());
    return;
}

void character_output::get_row(size_t r, Rcpp::StringVector::iterator out) { 
    get_row(r, out, 0, get_ncol());
    return;
}

void character_output::set_col(size_t c, Rcpp::StringVector::iterator out) { 
    set_col(c, out, 0, get_nrow());
    return;
}

void character_output::set_row(size_t r, Rcpp::StringVector::iterator out) { 
    set_row(r, out, 0, get_ncol());
    return;
}

/* Methods for the simple character matrix. */

simple_character_output::simple_character_output(size_t nr, size_t nc) : writer(nr, nc) {}

simple_character_output::~simple_character_output() {}

size_t simple_character_output::get_nrow() const {
    return writer.get_nrow();
}

size_t simple_character_output::get_ncol() const {
    return writer.get_ncol();
}

void simple_character_output::get_row(size_t r, Rcpp::StringVector::iterator out, size_t first, size_t last) { 
    writer.get_row(r, out, first, last);
    return;
}

void simple_character_output::get_col(size_t c, Rcpp::StringVector::iterator out, size_t first, size_t last) { 
    writer.get_col(c, out, first, last);
    return;
}

Rcpp::String simple_character_output::get(size_t r, size_t c) {
    return writer.get(r, c);
}

void simple_character_output::set_row(size_t r, Rcpp::StringVector::iterator in, size_t first, size_t last) { 
    writer.set_row(r, in, first, last);
    return;
}

void simple_character_output::set_col(size_t c, Rcpp::StringVector::iterator in, size_t first, size_t last) { 
    writer.set_col(c, in, first, last);
    return;
}

void simple_character_output::set(size_t r, size_t c, Rcpp::String in) {
    writer.set(r, c, in);
    return;
}

void simple_character_output::set_col_indexed(size_t c, size_t N, Rcpp::IntegerVector::iterator idx, Rcpp::StringVector::iterator val) {
    writer.set_col_indexed(c, N, idx, val);
    return;
}

void simple_character_output::set_row_indexed(size_t r, size_t N, Rcpp::IntegerVector::iterator idx, Rcpp::StringVector::iterator val) {
    writer.set_row_indexed(r, N, idx, val);
    return;
}

Rcpp::RObject simple_character_output::yield() {
    return writer.yield();
}

std::unique_ptr<character_output> simple_character_output::clone() const {
    return std::unique_ptr<character_output>(new simple_character_output(*this));
}

matrix_type simple_character_output::get_matrix_type() const {
    return writer.get_matrix_type();
}

/* Methods for the HDF5 output matrix. */

template<>
char HDF5_writer<char, STRSXP>::get_empty() { return '\0'; }

template<>
Rcpp::RObject HDF5_writer<char, STRSXP>::get_firstval() {
    std::vector<char> first(default_type.getSize());
    extract_one(0, 0, first.data());
    return Rcpp::StringVector::create(Rcpp::String(first.data()));
}

/* Methods for the HDF5 character matrix. */

HDF5_character_output::HDF5_character_output(size_t nr, size_t nc, size_t strlen, size_t chunk_nr, size_t chunk_nc, int compress) :
        bufsize(strlen+1), writer(nr, nc, chunk_nr, chunk_nc, compress, bufsize), 
        buffer(bufsize * std::max({ nr, nc, size_t(1) })) {}

HDF5_character_output::~HDF5_character_output() {}

size_t HDF5_character_output::get_nrow() const {
    return writer.get_nrow();
}

size_t HDF5_character_output::get_ncol() const {
    return writer.get_ncol();
}

void HDF5_character_output::get_row(size_t r, Rcpp::StringVector::iterator out, size_t first, size_t last) { 
    char* ref=buffer.data();
    writer.extract_row(r, ref, first, last);
    for (size_t c=first; c<last; ++c, ref+=bufsize, ++out) {
        (*out)=ref; 
    }
    return;
} 

void HDF5_character_output::get_col(size_t c, Rcpp::StringVector::iterator out, size_t first, size_t last) { 
    char* ref=buffer.data();
    writer.extract_col(c, ref, first, last);
    for (size_t r=first; r<last; ++r, ref+=bufsize, ++out) {
        (*out)=ref; 
    }
    return;
}
 
Rcpp::String HDF5_character_output::get(size_t r, size_t c) { 
    char* ref=buffer.data();
    writer.extract_one(r, c, ref);
    return ref;
}

void HDF5_character_output::set_row(size_t r, Rcpp::StringVector::iterator in, size_t first, size_t last) { 
    if (writer.get_ncol() + first >= last) { // ensure they can fit in 'buffer'; if not, it should trigger an error in insert_row().
        char* ref=buffer.data();
        for (size_t c=first; c<last; ++c, ref+=bufsize, ++in) {
            std::strncpy(ref, Rcpp::String(*in).get_cstring(), bufsize-1);
            ref[bufsize-1]='\0'; // strncpy only pads up to just before the last position.
        }
    }
    writer.insert_row(r, buffer.data(), first, last);
    return;
} 

void HDF5_character_output::set_col(size_t c, Rcpp::StringVector::iterator in, size_t first, size_t last) { 
    if (writer.get_nrow() + first >= last) { // ensure they can fit in 'buffer'.
        char* ref=buffer.data();
        for (size_t r=first; r<last; ++r, ref+=bufsize, ++in) {
            std::strncpy(ref, Rcpp::String(*in).get_cstring(), bufsize-1);
            ref[bufsize-1]='\0';
        }
    }
    writer.insert_col(c, buffer.data(), first, last);
    return;
}
 
void HDF5_character_output::set(size_t r, size_t c, Rcpp::String in) { 
    char* ref=buffer.data();
    std::strncpy(ref, in.get_cstring(), bufsize-1);
    ref[bufsize-1]='\0';
    writer.insert_one(r, c, ref);
    return;
}

void HDF5_character_output::set_col_indexed(size_t c, size_t N, Rcpp::IntegerVector::iterator idx, Rcpp::StringVector::iterator val) {
    if (buffer.size() < N*bufsize) {
        buffer.resize(N*bufsize);
    }

    char* ref=buffer.data();
    for (size_t i=0; i<N; ++i, ref+=bufsize, ++val) {
        std::strncpy(ref, Rcpp::String(*val).get_cstring(), bufsize-1);
        ref[bufsize-1]='\0';
    }
 
    writer.insert_col_indexed(c, N, idx, buffer.data());
    return;
}

void HDF5_character_output::set_row_indexed(size_t r, size_t N, Rcpp::IntegerVector::iterator idx, Rcpp::StringVector::iterator val) {
    if (buffer.size() < N*bufsize) {
        buffer.resize(N*bufsize);
    }

    char* ref=buffer.data();
    for (size_t i=0; i<N; ++i, ref+=bufsize, ++val) {
        std::strncpy(ref, Rcpp::String(*val).get_cstring(), bufsize-1);
        ref[bufsize-1]='\0';
    }
 
    writer.insert_row_indexed(r, N, idx, buffer.data());
    return;
}

Rcpp::RObject HDF5_character_output::yield() {
    return writer.yield();
}

std::unique_ptr<character_output> HDF5_character_output::clone() const {
    return std::unique_ptr<character_output>(new HDF5_character_output(*this));
}
 
matrix_type HDF5_character_output::get_matrix_type() const {
    return writer.get_matrix_type();
}

/* Dispatch definition */

std::unique_ptr<character_output> create_character_output(int nrow, int ncol, const output_param& param) {
    switch (param.get_mode()) {
        case SIMPLE:
            return std::unique_ptr<character_output>(new simple_character_output(nrow, ncol));
        case HDF5:
            return std::unique_ptr<character_output>(new HDF5_character_output(nrow, ncol,
                        param.get_strlen(), param.get_chunk_nrow(), param.get_chunk_ncol(), param.get_compression()));
        default:
            throw std::runtime_error("unsupported output mode for character matrices");
    }
}

}
