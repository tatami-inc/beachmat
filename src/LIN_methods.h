#ifndef BEACHMAT_LIN_METHODS_H
#define BEACHMAT_LIN_METHODS_H

namespace beachmat { 

/****************************************
 * Defining the common input interface. 
 ****************************************/

template<typename T, class V>
lin_matrix<T, V>::lin_matrix() {}

template<typename T, class V>
lin_matrix<T, V>::~lin_matrix() {}

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

template<typename T, class V, class M>
general_lin_matrix<T, V, M>::general_lin_matrix(const Rcpp::RObject& incoming) : mat(incoming) {}

template<typename T, class V, class M>
general_lin_matrix<T, V, M>::~general_lin_matrix() {}

template<typename T, class V, class M>
size_t general_lin_matrix<T, V, M>::get_nrow() const {
    return mat.get_nrow();
}

template<typename T, class V, class M>
size_t general_lin_matrix<T, V, M>::get_ncol() const {
    return mat.get_ncol();
}

template<typename T, class V, class M>
void general_lin_matrix<T, V, M>::get_col(size_t c, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
    mat.get_col(c, out, first, last);
    return;
}

template<typename T, class V, class M>
void general_lin_matrix<T, V, M>::get_col(size_t c, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
    mat.get_col(c, out, first, last);
    return;
}

template<typename T, class V, class M>
void general_lin_matrix<T, V, M>::get_row(size_t r, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
    mat.get_row(r, out, first, last);
    return;
}

template<typename T, class V, class M>
void general_lin_matrix<T, V, M>::get_row(size_t r, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
    mat.get_row(r, out, first, last);
    return;
}

template<typename T, class V, class M>
T general_lin_matrix<T, V, M>::get(size_t r, size_t c) {
    return mat.get(r, c);
}

template<typename T, class V, class M>
std::unique_ptr<lin_matrix<T, V> > general_lin_matrix<T, V, M>::clone() const {
    return std::unique_ptr<lin_matrix<T, V> >(new general_lin_matrix<T, V, M>(*this));
}

template<typename T, class V, class M> 
Rcpp::RObject general_lin_matrix<T, V, M>::yield() const {
    return mat.yield();
}

template<typename T, class V, class M> 
matrix_type general_lin_matrix<T, V, M>::get_matrix_type() const {
    return mat.get_matrix_type();
}

/* Defining specific interface for simple matrices. */

template <typename T, class V>
simple_lin_matrix<T, V>::simple_lin_matrix(const Rcpp::RObject& in) : simple_lin_precursor<T, V>(in) {}

template <typename T, class V>
simple_lin_matrix<T, V>::~simple_lin_matrix() {} 

template <typename T, class V>
typename V::iterator simple_lin_matrix<T, V>::get_const_col(size_t c, typename V::iterator work, size_t first, size_t last) {
    return this->mat.get_const_col(c, first, last);
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
    return this->mat.get_const_col(c, first, last);
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
    size_t nzero=this->mat.get_const_col_nonzero(c, iIt, out, first, last);
    return const_col_indexed_info<V>(nzero, iIt, out); 
}

template <typename T, class V>
std::unique_ptr<lin_matrix<T, V> > Csparse_lin_matrix<T, V>::clone() const {
    return std::unique_ptr<lin_matrix<T, V> >(new Csparse_lin_matrix<T, V>(*this));
}

/* Defining the HDF5 interface (or specifically, the base object for LIN data). */

template<typename T, int RTYPE>
HDF5_lin_base<T, RTYPE>::HDF5_lin_base(const Rcpp::RObject& incoming) : HDF5_matrix<T, RTYPE>(incoming) {}

template<typename T, int RTYPE>
HDF5_lin_base<T, RTYPE>::~HDF5_lin_base() {}

template<typename T, int RTYPE>
void HDF5_lin_base<T, RTYPE>::get_col(size_t c, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
    this->extract_col(c, &(*out), H5::PredType::NATIVE_INT32, first, last);
    return;
}

template<typename T, int RTYPE>
void HDF5_lin_base<T, RTYPE>::get_col(size_t c, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
    this->extract_col(c, &(*out), H5::PredType::NATIVE_DOUBLE, first, last);
    return;
}

template<typename T, int RTYPE>
void HDF5_lin_base<T, RTYPE>::get_row(size_t r, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
    this->extract_row(r, &(*out), H5::PredType::NATIVE_INT32, first, last);
    return;
}

template<typename T, int RTYPE>
void HDF5_lin_base<T, RTYPE>::get_row(size_t r, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
    this->extract_row(r, &(*out), H5::PredType::NATIVE_DOUBLE, first, last);
    return;
}

/* Defining the DelayedMatrix interface. */

template<typename T, class V>
delayed_lin_matrix<T, V>::delayed_lin_matrix(const Rcpp::RObject& incoming) : original(incoming), seed_ptr(nullptr) {
    check_DelayedMatrix(incoming);

    // Trying to generate the seed, if it's a valid object in itself.
    if (only_delayed_coord_changes(incoming)) {
        seed_ptr=generate_seed(incoming);
    }
        
    // If the seed is still NULL, we switch to a chunked matrix format.
    if (seed_ptr.get()==NULL) { 
        seed_ptr=std::unique_ptr<lin_matrix<T, V> >(new delayed_lin_matrix<T, V>::enslaved(incoming));
        transformer=delayed_coord_transformer<T, V>(seed_ptr.get());
    } else {
        transformer=delayed_coord_transformer<T, V>(incoming, seed_ptr.get());
    }

    return;
}

template<typename T, class V> 
delayed_lin_matrix<T, V>::~delayed_lin_matrix() {}

template<typename T, class V> 
delayed_lin_matrix<T, V>::delayed_lin_matrix(const delayed_lin_matrix<T, V>& other) : original(other.original), 
        seed_ptr(other.seed_ptr->clone()), transformer(other.transformer) {}

template<typename T, class V> 
delayed_lin_matrix<T, V>& delayed_lin_matrix<T, V>::operator=(const delayed_lin_matrix<T, V>& other) {
    original=other.original;
    seed_ptr=other.seed_ptr->clone();
    transformer=other.transformer;
    return *this;
}

template<typename T, class V>
size_t delayed_lin_matrix<T, V>::get_nrow() const {
    return transformer.get_nrow();
}

template<typename T, class V>
size_t delayed_lin_matrix<T, V>::get_ncol() const {
    return transformer.get_ncol();
}

template<typename T, class V>
void delayed_lin_matrix<T, V>::get_col(size_t c, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
    transformer.get_col(seed_ptr.get(), c, out, first, last);
    return;
}

template<typename T, class V>
void delayed_lin_matrix<T, V>::get_col(size_t c, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
    transformer.get_col(seed_ptr.get(), c, out, first, last);
    return;
}

template<typename T, class V>
void delayed_lin_matrix<T, V>::get_row(size_t r, Rcpp::IntegerVector::iterator out, size_t first, size_t last) {
    transformer.get_row(seed_ptr.get(), r, out, first, last);
    return;
}

template<typename T, class V>
void delayed_lin_matrix<T, V>::get_row(size_t r, Rcpp::NumericVector::iterator out, size_t first, size_t last) {
    transformer.get_row(seed_ptr.get(), r, out, first, last);
    return;
}

template<typename T, class V>
T delayed_lin_matrix<T, V>::get(size_t r, size_t c) {
    return transformer.get(seed_ptr.get(), r, c);
}

template<typename T, class V>
std::unique_ptr<lin_matrix<T, V> > delayed_lin_matrix<T, V>::clone() const {
    return std::unique_ptr<lin_matrix<T, V> >(new delayed_lin_matrix<T, V>(*this));
}

template<typename T, class V> 
Rcpp::RObject delayed_lin_matrix<T, V>::yield() const {
    return original;
}

template<typename T, class V>
matrix_type delayed_lin_matrix<T, V>::get_matrix_type() const { 
    return DELAYED;
}

}

#endif
