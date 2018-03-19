#ifndef BEACHMAT_LIN_OUTFUN_H
#define BEACHMAT_LIN_OUTFUN_H

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

template<typename T, class V>
typename V::iterator lin_output<T, V>::get_const_col(size_t c, typename V::iterator work) {
    return get_const_col(c, work, 0, get_nrow());
}

template<typename T, class V>
typename V::iterator lin_output<T, V>::get_const_col(size_t c, typename V::iterator work, size_t first, size_t last) {
    get_col(c, work, first, last);
    return work;
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

template<typename T, class V, class M>
general_lin_output<T, V, M>::general_lin_output(size_t nr, size_t nc) : mat(nr, nc) {}

template<typename T, class V, class M>
general_lin_output<T, V, M>::~general_lin_output() {}

// Getters:
template<typename T, class V, class M>
size_t general_lin_output<T, V, M>::get_nrow() const {
    return mat.get_nrow();
}

template<typename T, class V, class M>
size_t general_lin_output<T, V, M>::get_ncol() const {
    return mat.get_ncol();
}

template<typename T, class V, class M>
void general_lin_output<T, V, M>::get_col(size_t c, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
    mat.get_col(c, out, first, last);
    return;
}

template<typename T, class V, class M>
void general_lin_output<T, V, M>::get_col(size_t c, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
    mat.get_col(c, out, first, last);
    return;
}

template<typename T, class V, class M>
void general_lin_output<T, V, M>::get_row(size_t r, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
    mat.get_row(r, out, first, last);
    return;
}

template<typename T, class V, class M>
void general_lin_output<T, V, M>::get_row(size_t r, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
    mat.get_row(r, out, first, last);
    return;
}

template<typename T, class V, class M>
T general_lin_output<T, V, M>::get(size_t r, size_t c) {
    return mat.get(r, c);
}

// Setters:
template<typename T, class V, class M>
void general_lin_output<T, V, M>::set_col(size_t c, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
    mat.set_col(c, out, first, last);
    return;
}

template<typename T, class V, class M>
void general_lin_output<T, V, M>::set_col(size_t c, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
    mat.set_col(c, out, first, last);
    return;
}

template<typename T, class V, class M>
void general_lin_output<T, V, M>::set_row(size_t r, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
    mat.set_row(r, out, first, last);
    return;
}

template<typename T, class V, class M>
void general_lin_output<T, V, M>::set_row(size_t r, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
    mat.set_row(r, out, first, last);
    return;
}

template<typename T, class V, class M>
void general_lin_output<T, V, M>::set(size_t r, size_t c, T in) {
    mat.set(r, c, in);
    return;
}

template<typename T, class V, class M>
void general_lin_output<T, V, M>::set_col_indexed(size_t c, const const_col_indexed_info<Rcpp::IntegerVector>& info) {
    mat.set_col_indexed(c, std::get<0>(info), std::get<1>(info), std::get<2>(info));
    return; 
}

template<typename T, class V, class M>
void general_lin_output<T, V, M>::set_col_indexed(size_t c, const const_col_indexed_info<Rcpp::NumericVector>& info) {
    mat.set_col_indexed(c, std::get<0>(info), std::get<1>(info), std::get<2>(info));
    return; 
}

// Other functions:
template<typename T, class V, class M>
Rcpp::RObject general_lin_output<T, V, M>::yield() {
    return mat.yield();
}

template<typename T, class V, class M>
std::unique_ptr<lin_output<T, V> > general_lin_output<T, V, M>::clone() const {
    return std::unique_ptr<lin_output<T, V> >(new general_lin_output<T, V, M>(*this));
}

template<typename T, class V, class M>
matrix_type general_lin_output<T, V, M>::get_matrix_type() const {
    return mat.get_matrix_type();
}

/* Defining the simple output interface. */ 

template<typename T, class V>
simple_lin_output<T, V>::simple_lin_output(size_t nr, size_t nc) : simple_lin_output_precursor<T, V>(nr, nc) {}

template<typename T, class V>
simple_lin_output<T, V>::~simple_lin_output() {}

template<typename T, class V>
typename V::iterator simple_lin_output<T, V>::get_const_col(size_t c, typename V::iterator work, size_t first, size_t last) {
    return this->mat.get_const_col(c, first, last);
}

template<typename T, class V>
std::unique_ptr<lin_output<T, V> > simple_lin_output<T, V>::clone() const {
    return std::unique_ptr<lin_output<T, V> >(new simple_lin_output<T, V>(*this));
}

/* Defining the sparse output interface. */ 

template<typename T, class V>
sparse_lin_output<T, V>::sparse_lin_output(size_t nr, size_t nc) : sparse_lin_output_precursor<T, V>(nr, nc) {}

template<typename T, class V>
sparse_lin_output<T, V>::~sparse_lin_output() {}

template<typename T, class V>
std::unique_ptr<lin_output<T, V> > sparse_lin_output<T, V>::clone() const {
    return std::unique_ptr<lin_output<T, V> >(new sparse_lin_output<T, V>(*this));
}

/* Defining the HDF5 output interface. */

template<typename T, class V>
HDF5_lin_output<T, V>::HDF5_lin_output(size_t nr, size_t nc, size_t chunk_nr, size_t chunk_nc, int compress) : 
    mat(nr, nc, chunk_nr, chunk_nc, compress) {}

template<typename T, class V>
HDF5_lin_output<T, V>::~HDF5_lin_output() {}

template<typename T, class V>
size_t HDF5_lin_output<T, V>::get_nrow() const {
    return mat.get_nrow();
}

template<typename T, class V>
size_t HDF5_lin_output<T, V>::get_ncol() const {
    return mat.get_ncol();
}

template<typename T, class V>
T HDF5_lin_output<T, V>::get(size_t r, size_t c) {
    T out;
    mat.extract_one(r, c, &out);
    return out;
}

template<typename T, class V>
void HDF5_lin_output<T, V>::set(size_t r, size_t c, T in) {
    mat.insert_one(r, c, &in);
    return;
}

template<typename T, class V>
void HDF5_lin_output<T, V>::get_row(size_t r, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
    mat.extract_row(r, &(*out), H5::PredType::NATIVE_INT32, first, last);
    return;
}

template<typename T, class V>
void HDF5_lin_output<T, V>::get_row(size_t r, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
    mat.extract_row(r, &(*out), H5::PredType::NATIVE_DOUBLE, first, last);
    return;
}

template<typename T, class V>
void HDF5_lin_output<T, V>::get_col(size_t c, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
    mat.extract_col(c, &(*out), H5::PredType::NATIVE_INT32, first, last);
    return;
}

template<typename T, class V>
void HDF5_lin_output<T, V>::get_col(size_t c, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
    mat.extract_col(c, &(*out), H5::PredType::NATIVE_DOUBLE, first, last);
    return;
}

template<typename T, class V>
void HDF5_lin_output<T, V>::set_row(size_t r, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
    mat.insert_row(r, &(*out), H5::PredType::NATIVE_INT32, first, last);
    return;
}

template<typename T, class V>
void HDF5_lin_output<T, V>::set_row(size_t r, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
    mat.insert_row(r, &(*out), H5::PredType::NATIVE_DOUBLE, first, last);
    return;
}

template<typename T, class V>
void HDF5_lin_output<T, V>::set_col(size_t c, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
    mat.insert_col(c, &(*out), H5::PredType::NATIVE_INT32, first, last);
    return;
}

template<typename T, class V>
void HDF5_lin_output<T, V>::set_col(size_t c, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
    mat.insert_col(c, &(*out), H5::PredType::NATIVE_DOUBLE, first, last);
    return;
}

template<typename T, class V>
void HDF5_lin_output<T, V>::set_col_indexed(size_t c, const const_col_indexed_info<Rcpp::IntegerVector>& info) {
    size_t N=std::get<0>(info);
    auto idx=std::get<1>(info);
    auto val=std::get<2>(info);
    for (size_t i=0; i<N; ++i, ++idx, ++val) {
        T tmp=*val;
        mat.insert_one(*idx, c, &tmp);
    }
    return;
}
 
template<typename T, class V>
void HDF5_lin_output<T, V>::set_col_indexed(size_t c, const const_col_indexed_info<Rcpp::NumericVector>& info) {
    size_t N=std::get<0>(info);
    auto idx=std::get<1>(info);
    auto val=std::get<2>(info);
    for (size_t i=0; i<N; ++i, ++idx, ++val) {
        T tmp=*val;
        mat.insert_one(*idx, c, &tmp);
    }
    return;
}

template<typename T, class V>
Rcpp::RObject HDF5_lin_output<T, V>::yield() {
    return mat.yield();
}

template<typename T, class V>
std::unique_ptr<lin_output<T, V> > HDF5_lin_output<T, V>::clone() const {
    return std::unique_ptr<lin_output<T, V> >(new HDF5_lin_output<T, V>(*this));
}

template<typename T, class V>
matrix_type HDF5_lin_output<T, V>::get_matrix_type() const {
    return mat.get_matrix_type();
}

}

#endif
