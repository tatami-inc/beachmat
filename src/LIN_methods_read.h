#ifndef BEACHMAT_LIN_METHODS_READ_H
#define BEACHMAT_LIN_METHODS_READ_H

namespace beachmat { 

/****************************************
 * Defining the common input interface. 
 ****************************************/

template<typename T, class V>
lin_matrix<T, V>::lin_matrix() {}

template<typename T, class V>
lin_matrix<T, V>::~lin_matrix() {}

// Basic getters.

template<typename T, class V>
void lin_matrix<T, V>::get_col(size_t c, Rcpp::IntegerVector::iterator out) {
    get_col(c, out, 0, get_nrow());
    return;
}

template<typename T, class V>
void lin_matrix<T, V>::get_col(size_t c, Rcpp::NumericVector::iterator out) {
    get_col(c, out, 0, get_nrow());
    return;
}

template<typename T, class V>
void lin_matrix<T, V>::get_row(size_t r, Rcpp::IntegerVector::iterator out) {
    get_row(r, out, 0, get_ncol());
    return;
}

template<typename T, class V>
void lin_matrix<T, V>::get_row(size_t r, Rcpp::NumericVector::iterator out) {
    get_row(r, out, 0, get_ncol());
    return;
}

// Multi getters.

template<typename T, class V>
void lin_matrix<T, V>::get_cols(Rcpp::IntegerVector::iterator it, size_t n, Rcpp::IntegerVector::iterator out) {
    get_cols(it, n, out, 0, get_nrow());
    return;
}

template<typename T, class V>
void lin_matrix<T, V>::get_cols(Rcpp::IntegerVector::iterator it, size_t n, Rcpp::NumericVector::iterator out) {
    get_cols(it, n, out, 0, get_nrow());
    return;
}

template<typename T, class V>
void lin_matrix<T, V>::get_rows(Rcpp::IntegerVector::iterator it, size_t n, Rcpp::IntegerVector::iterator out) {
    get_rows(it, n, out, 0, get_ncol());
    return;
}

template<typename T, class V>
void lin_matrix<T, V>::get_rows(Rcpp::IntegerVector::iterator it, size_t n, Rcpp::NumericVector::iterator out) {
    get_rows(it, n, out, 0, get_ncol());
    return;
}

// Specialist getters.

template<typename T, class V>
typename V::iterator lin_matrix<T, V>::get_const_col(size_t c, typename V::iterator work) {
    return get_const_col(c, work, 0, get_nrow());
}

template<typename T, class V>
typename V::iterator lin_matrix<T, V>::get_const_col(size_t c, typename V::iterator work, size_t first, size_t last) {
    get_col(c, work, first, last);
    return work;
}

template<typename T, class V>
const_col_indexed_info<V> lin_matrix<T, V>::get_const_col_indexed(size_t c, typename V::iterator work) {
    return get_const_col_indexed(c, work, 0, get_nrow());
}

template<typename T, class V>
const_col_indexed_info<V> lin_matrix<T, V>::get_const_col_indexed(size_t c, typename V::iterator work, size_t first, size_t last) {
    if (indices.size()!=this->get_nrow()) {
        indices=Rcpp::IntegerVector(this->get_nrow());
        std::iota(indices.begin(), indices.end(), 0); // populating with indices.
    }
    return const_col_indexed_info<V>(last - first, indices.begin() + first, get_const_col(c, work, first, last));
}

/* Defining the general interface. */

template<typename T, class V, class RDR>
general_lin_matrix<T, V, RDR>::general_lin_matrix(const Rcpp::RObject& incoming) : reader(incoming) {}

template<typename T, class V, class RDR>
general_lin_matrix<T, V, RDR>::~general_lin_matrix() {}

template<typename T, class V, class RDR>
size_t general_lin_matrix<T, V, RDR>::get_nrow() const {
    return reader.get_nrow();
}

template<typename T, class V, class RDR>
size_t general_lin_matrix<T, V, RDR>::get_ncol() const {
    return reader.get_ncol();
}

// Basic getters.

template<typename T, class V, class RDR>
void general_lin_matrix<T, V, RDR>::get_col(size_t c, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
    reader.get_col(c, out, first, last);
    return;
}

template<typename T, class V, class RDR>
void general_lin_matrix<T, V, RDR>::get_col(size_t c, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
    reader.get_col(c, out, first, last);
    return;
}

template<typename T, class V, class RDR>
void general_lin_matrix<T, V, RDR>::get_row(size_t r, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
    reader.get_row(r, out, first, last);
    return;
}

template<typename T, class V, class RDR>
void general_lin_matrix<T, V, RDR>::get_row(size_t r, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
    reader.get_row(r, out, first, last);
    return;
}

template<typename T, class V, class RDR>
T general_lin_matrix<T, V, RDR>::get(size_t r, size_t c) {
    return reader.get(r, c);
}

// Multi getters.

template<typename T, class V, class RDR>
void general_lin_matrix<T, V, RDR>::get_cols(Rcpp::IntegerVector::iterator it, size_t n, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
    reader.get_cols(it, n, out, first, last);
    return;
}

template<typename T, class V, class RDR>
void general_lin_matrix<T, V, RDR>::get_cols(Rcpp::IntegerVector::iterator it, size_t n, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
    reader.get_cols(it, n, out, first, last);
    return;
}

template<typename T, class V, class RDR>
void general_lin_matrix<T, V, RDR>::get_rows(Rcpp::IntegerVector::iterator it, size_t n, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
    reader.get_rows(it, n, out, first, last);
    return;
}

template<typename T, class V, class RDR>
void general_lin_matrix<T, V, RDR>::get_rows(Rcpp::IntegerVector::iterator it, size_t n, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
    reader.get_rows(it, n, out, first, last);
    return;
}

// Other methods.

template<typename T, class V, class RDR>
std::unique_ptr<lin_matrix<T, V> > general_lin_matrix<T, V, RDR>::clone() const {
    return std::unique_ptr<lin_matrix<T, V> >(new general_lin_matrix<T, V, RDR>(*this));
}

template<typename T, class V, class RDR> 
Rcpp::RObject general_lin_matrix<T, V, RDR>::yield() const {
    return reader.yield();
}

template<typename T, class V, class RDR> 
matrix_type general_lin_matrix<T, V, RDR>::get_matrix_type() const {
    return reader.get_matrix_type();
}

/* Defining specific interface for simple matrices. */

template <typename T, class V>
simple_lin_matrix<T, V>::simple_lin_matrix(const Rcpp::RObject& in) : simple_lin_precursor<T, V>(in) {}

template <typename T, class V>
simple_lin_matrix<T, V>::~simple_lin_matrix() {} 

template <typename T, class V>
typename V::iterator simple_lin_matrix<T, V>::get_const_col(size_t c, typename V::iterator work, size_t first, size_t last) {
    return this->reader.get_const_col(c, first, last);
}

template <typename T, class V>
std::unique_ptr<lin_matrix<T, V> > simple_lin_matrix<T, V>::clone() const {
    return std::unique_ptr<lin_matrix<T, V> >(new simple_lin_matrix<T, V>(*this));
}

/* Defining specific interface for dense matrices. */

template <typename T, class V>
dense_lin_matrix<T, V>::dense_lin_matrix(const Rcpp::RObject& in) : dense_lin_precursor<T, V>(in) {}

template <typename T, class V>
dense_lin_matrix<T, V>::~dense_lin_matrix() {} 

template <typename T, class V>
typename V::iterator dense_lin_matrix<T, V>::get_const_col(size_t c, typename V::iterator work, size_t first, size_t last) {
    return this->reader.get_const_col(c, first, last);
}

template <typename T, class V>
std::unique_ptr<lin_matrix<T, V> > dense_lin_matrix<T, V>::clone() const {
    return std::unique_ptr<lin_matrix<T, V> >(new dense_lin_matrix<T, V>(*this));
}

/* Defining specific interface for sparse matrices. */

template <typename T, class V>
Csparse_lin_matrix<T, V>::Csparse_lin_matrix(const Rcpp::RObject& in) : Csparse_lin_precursor<T, V>(in) {}

template <typename T, class V>
Csparse_lin_matrix<T, V>::~Csparse_lin_matrix() {} 

template <typename T, class V>
const_col_indexed_info<V> Csparse_lin_matrix<T, V>::get_const_col_indexed(size_t c, typename V::iterator out, size_t first, size_t last) {
    Rcpp::IntegerVector::iterator iIt;
    size_t nzero=this->reader.get_const_col_nonzero(c, iIt, out, first, last);
    return const_col_indexed_info<V>(nzero, iIt, out); 
}

template <typename T, class V>
std::unique_ptr<lin_matrix<T, V> > Csparse_lin_matrix<T, V>::clone() const {
    return std::unique_ptr<lin_matrix<T, V> >(new Csparse_lin_matrix<T, V>(*this));
}

/* Defining the helper class contained inside the HDF5 interface. */

template<typename T, int RTYPE>
HDF5_lin_reader<T, RTYPE>::HDF5_lin_reader(const Rcpp::RObject& incoming) : HDF5_reader<T, RTYPE>(incoming) {}

template<typename T, int RTYPE>
HDF5_lin_reader<T, RTYPE>::~HDF5_lin_reader() {}

template<typename T, int RTYPE>
void HDF5_lin_reader<T, RTYPE>::get_col(size_t c, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
    this->extract_col(c, static_cast<int*>(out), H5::PredType::NATIVE_INT32, first, last);
    return;
}

template<typename T, int RTYPE>
void HDF5_lin_reader<T, RTYPE>::get_col(size_t c, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
    this->extract_col(c, static_cast<double*>(out), H5::PredType::NATIVE_DOUBLE, first, last);
    return;
}

template<typename T, int RTYPE>
void HDF5_lin_reader<T, RTYPE>::get_row(size_t r, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
    this->extract_row(r, static_cast<int*>(out), H5::PredType::NATIVE_INT32, first, last);
    return;
}

template<typename T, int RTYPE>
void HDF5_lin_reader<T, RTYPE>::get_row(size_t r, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
    this->extract_row(r, static_cast<double*>(out), H5::PredType::NATIVE_DOUBLE, first, last);
    return;
}

// Multi getters.

template<typename T, int RTYPE>
void HDF5_lin_reader<T, RTYPE>::get_cols(Rcpp::IntegerVector::iterator it, size_t n, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
    this->extract_cols(it, n, static_cast<int*>(out), H5::PredType::NATIVE_INT32, first, last);
    return;
}

template<typename T, int RTYPE>
void HDF5_lin_reader<T, RTYPE>::get_cols(Rcpp::IntegerVector::iterator it, size_t n, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
    this->extract_cols(it, n, static_cast<double*>(out), H5::PredType::NATIVE_DOUBLE, first, last);
    return;
}

template<typename T, int RTYPE>
void HDF5_lin_reader<T, RTYPE>::get_rows(Rcpp::IntegerVector::iterator it, size_t n, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
    this->extract_rows(it, n, static_cast<int*>(out), H5::PredType::NATIVE_INT32, first, last);
    return;
}

template<typename T, int RTYPE>
void HDF5_lin_reader<T, RTYPE>::get_rows(Rcpp::IntegerVector::iterator it, size_t n, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
    this->extract_rows(it, n, static_cast<double*>(out), H5::PredType::NATIVE_DOUBLE, first, last);
    return;
}

}

#endif
