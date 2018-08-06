#ifndef BEACHMAT_LIN_METHODS_WRITE_H
#define BEACHMAT_LIN_METHODS_WRITE_H

namespace beachmat { 

/****************************************
 * Defining the common output interface. 
 ****************************************/

template<typename T, class V>
lin_output<T, V>::lin_output() {}

template<typename T, class V>
lin_output<T, V>::~lin_output() {}

// Getters:
template<typename T, class V>
void lin_output<T, V>::get_col(size_t c, Rcpp::IntegerVector::iterator out) {
    get_col(c, out, 0, get_nrow());
    return;
}

template<typename T, class V>
void lin_output<T, V>::get_col(size_t c, Rcpp::NumericVector::iterator out) {
    get_col(c, out, 0, get_nrow());
    return;
}

template<typename T, class V>
void lin_output<T, V>::get_row(size_t r, Rcpp::IntegerVector::iterator out) {
    get_row(r, out, 0, get_ncol());
    return;
}

template<typename T, class V>
void lin_output<T, V>::get_row(size_t r, Rcpp::NumericVector::iterator out) {
    get_row(r, out, 0, get_ncol());
    return;
}

// Setters:
template<typename T, class V>
void lin_output<T, V>::set_col(size_t c, Rcpp::IntegerVector::iterator out) {
    set_col(c, out, 0, get_nrow());
    return;
}

template<typename T, class V>
void lin_output<T, V>::set_col(size_t c, Rcpp::NumericVector::iterator out) {
    set_col(c, out, 0, get_nrow());
    return;
}

template<typename T, class V>
void lin_output<T, V>::set_row(size_t r, Rcpp::IntegerVector::iterator out) {
    set_row(r, out, 0, get_ncol());
    return;
}

template<typename T, class V>
void lin_output<T, V>::set_row(size_t r, Rcpp::NumericVector::iterator out) {
    set_row(r, out, 0, get_ncol());
    return;
}

/* Defining the general output interface. */ 

template<typename T, class V, class WTR>
general_lin_output<T, V, WTR>::general_lin_output(size_t nr, size_t nc) : writer(nr, nc) {}

template<typename T, class V, class WTR>
general_lin_output<T, V, WTR>::~general_lin_output() {}

// Getters:
template<typename T, class V, class WTR>
size_t general_lin_output<T, V, WTR>::get_nrow() const {
    return writer.get_nrow();
}

template<typename T, class V, class WTR>
size_t general_lin_output<T, V, WTR>::get_ncol() const {
    return writer.get_ncol();
}

template<typename T, class V, class WTR>
void general_lin_output<T, V, WTR>::get_col(size_t c, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
    writer.get_col(c, out, first, last);
    return;
}

template<typename T, class V, class WTR>
void general_lin_output<T, V, WTR>::get_col(size_t c, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
    writer.get_col(c, out, first, last);
    return;
}

template<typename T, class V, class WTR>
void general_lin_output<T, V, WTR>::get_row(size_t r, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
    writer.get_row(r, out, first, last);
    return;
}

template<typename T, class V, class WTR>
void general_lin_output<T, V, WTR>::get_row(size_t r, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
    writer.get_row(r, out, first, last);
    return;
}

template<typename T, class V, class WTR>
T general_lin_output<T, V, WTR>::get(size_t r, size_t c) {
    return writer.get(r, c);
}

// Setters:
template<typename T, class V, class WTR>
void general_lin_output<T, V, WTR>::set_col(size_t c, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
    writer.set_col(c, out, first, last);
    return;
}

template<typename T, class V, class WTR>
void general_lin_output<T, V, WTR>::set_col(size_t c, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
    writer.set_col(c, out, first, last);
    return;
}

template<typename T, class V, class WTR>
void general_lin_output<T, V, WTR>::set_row(size_t r, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
    writer.set_row(r, out, first, last);
    return;
}

template<typename T, class V, class WTR>
void general_lin_output<T, V, WTR>::set_row(size_t r, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
    writer.set_row(r, out, first, last);
    return;
}

template<typename T, class V, class WTR>
void general_lin_output<T, V, WTR>::set(size_t r, size_t c, T in) {
    writer.set(r, c, in);
    return;
}

template<typename T, class V, class WTR>
void general_lin_output<T, V, WTR>::set_col_indexed(size_t c, size_t N, Rcpp::IntegerVector::iterator idx, Rcpp::IntegerVector::iterator val) {
    writer.set_col_indexed(c, N, idx, val);
    return; 
}

template<typename T, class V, class WTR>
void general_lin_output<T, V, WTR>::set_col_indexed(size_t c, size_t N, Rcpp::IntegerVector::iterator idx, Rcpp::NumericVector::iterator val) {
    writer.set_col_indexed(c, N, idx, val);
    return; 
}

template<typename T, class V, class WTR>
void general_lin_output<T, V, WTR>::set_row_indexed(size_t r, size_t N, Rcpp::IntegerVector::iterator idx, Rcpp::IntegerVector::iterator val) {
    writer.set_row_indexed(r, N, idx, val);
    return; 
}

template<typename T, class V, class WTR>
void general_lin_output<T, V, WTR>::set_row_indexed(size_t r, size_t N, Rcpp::IntegerVector::iterator idx, Rcpp::NumericVector::iterator val) {
    writer.set_row_indexed(r, N, idx, val);
    return; 
}

// Other functions:
template<typename T, class V, class WTR>
Rcpp::RObject general_lin_output<T, V, WTR>::yield() {
    return writer.yield();
}

template<typename T, class V, class WTR>
std::unique_ptr<lin_output<T, V> > general_lin_output<T, V, WTR>::clone() const {
    return std::unique_ptr<lin_output<T, V> >(new general_lin_output<T, V, WTR>(*this));
}

template<typename T, class V, class WTR>
matrix_type general_lin_output<T, V, WTR>::get_matrix_type() const {
    return writer.get_matrix_type();
}

/* Defining the simple output interface. */ 

template<typename T, class V>
simple_lin_output<T, V>::simple_lin_output(size_t nr, size_t nc) : simple_lin_output_precursor<T, V>(nr, nc) {}

template<typename T, class V>
simple_lin_output<T, V>::~simple_lin_output() {}

template<typename T, class V>
typename V::iterator simple_lin_output<T, V>::get_const_col(size_t c, typename V::iterator work, size_t first, size_t last) {
    return this->writer.get_const_col(c, first, last);
}

template<typename T, class V>
std::unique_ptr<lin_output<T, V> > simple_lin_output<T, V>::clone() const {
    return std::unique_ptr<lin_output<T, V> >(new simple_lin_output<T, V>(*this));
}

/* Defining the HDF5 output interface. */

template<typename T, class V, int RTYPE>
HDF5_lin_output<T, V, RTYPE>::HDF5_lin_output(size_t nr, size_t nc, size_t chunk_nr, size_t chunk_nc, int compress) : 
    writer(nr, nc, chunk_nr, chunk_nc, compress) {}

template<typename T, class V, int RTYPE>
HDF5_lin_output<T, V, RTYPE>::~HDF5_lin_output() {}

template<typename T, class V, int RTYPE>
size_t HDF5_lin_output<T, V, RTYPE>::get_nrow() const {
    return writer.get_nrow();
}

template<typename T, class V, int RTYPE>
size_t HDF5_lin_output<T, V, RTYPE>::get_ncol() const {
    return writer.get_ncol();
}

template<typename T, class V, int RTYPE>
T HDF5_lin_output<T, V, RTYPE>::get(size_t r, size_t c) {
    T out;
    writer.extract_one(r, c, &out);
    return out;
}

template<typename T, class V, int RTYPE>
void HDF5_lin_output<T, V, RTYPE>::set(size_t r, size_t c, T in) {
    writer.insert_one(r, c, &in);
    return;
}

template<typename T, class V, int RTYPE>
void HDF5_lin_output<T, V, RTYPE>::get_row(size_t r, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
    writer.extract_row(r, &(*out), H5::PredType::NATIVE_INT32, first, last);
    return;
}

template<typename T, class V, int RTYPE>
void HDF5_lin_output<T, V, RTYPE>::get_row(size_t r, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
    writer.extract_row(r, &(*out), H5::PredType::NATIVE_DOUBLE, first, last);
    return;
}

template<typename T, class V, int RTYPE>
void HDF5_lin_output<T, V, RTYPE>::get_col(size_t c, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
    writer.extract_col(c, &(*out), H5::PredType::NATIVE_INT32, first, last);
    return;
}

template<typename T, class V, int RTYPE>
void HDF5_lin_output<T, V, RTYPE>::get_col(size_t c, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
    writer.extract_col(c, &(*out), H5::PredType::NATIVE_DOUBLE, first, last);
    return;
}

template<typename T, class V, int RTYPE>
void HDF5_lin_output<T, V, RTYPE>::set_row(size_t r, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
    writer.insert_row(r, &(*out), H5::PredType::NATIVE_INT32, first, last);
    return;
}

template<typename T, class V, int RTYPE>
void HDF5_lin_output<T, V, RTYPE>::set_row(size_t r, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
    writer.insert_row(r, &(*out), H5::PredType::NATIVE_DOUBLE, first, last);
    return;
}

template<typename T, class V, int RTYPE>
void HDF5_lin_output<T, V, RTYPE>::set_col(size_t c, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
    writer.insert_col(c, &(*out), H5::PredType::NATIVE_INT32, first, last);
    return;
}

template<typename T, class V, int RTYPE>
void HDF5_lin_output<T, V, RTYPE>::set_col(size_t c, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
    writer.insert_col(c, &(*out), H5::PredType::NATIVE_DOUBLE, first, last);
    return;
}

template<typename T, class V, int RTYPE>
void HDF5_lin_output<T, V, RTYPE>::set_col_indexed(size_t c, size_t N, Rcpp::IntegerVector::iterator idx, Rcpp::IntegerVector::iterator val) {
    writer.insert_col_indexed(c, N, idx, val, H5::PredType::NATIVE_INT32);
    return;
}
 
template<typename T, class V, int RTYPE>
void HDF5_lin_output<T, V, RTYPE>::set_col_indexed(size_t c, size_t N, Rcpp::IntegerVector::iterator idx, Rcpp::NumericVector::iterator val) {
    writer.insert_col_indexed(c, N, idx, val, H5::PredType::NATIVE_DOUBLE);
    return;
}

template<typename T, class V, int RTYPE>
void HDF5_lin_output<T, V, RTYPE>::set_row_indexed(size_t r, size_t N, Rcpp::IntegerVector::iterator idx, Rcpp::IntegerVector::iterator val) {
    writer.insert_row_indexed(r, N, idx, val, H5::PredType::NATIVE_INT32);
    return;
}
 
template<typename T, class V, int RTYPE>
void HDF5_lin_output<T, V, RTYPE>::set_row_indexed(size_t r, size_t N, Rcpp::IntegerVector::iterator idx, Rcpp::NumericVector::iterator val) {
    writer.insert_row_indexed(r, N, idx, val, H5::PredType::NATIVE_DOUBLE);
    return;
}

template<typename T, class V, int RTYPE>
Rcpp::RObject HDF5_lin_output<T, V, RTYPE>::yield() {
    return writer.yield();
}

template<typename T, class V, int RTYPE>
std::unique_ptr<lin_output<T, V> > HDF5_lin_output<T, V, RTYPE>::clone() const {
    return std::unique_ptr<lin_output<T, V> >(new HDF5_lin_output<T, V, RTYPE>(*this));
}

template<typename T, class V, int RTYPE>
matrix_type HDF5_lin_output<T, V, RTYPE>::get_matrix_type() const {
    return writer.get_matrix_type();
}

}

#endif
