#ifndef UNKNOWN_MATRIX_H
#define UNKNOWN_MATRIX_H

namespace beachmat {

/* The 'unknown_reader' class will realize chunks of the input RObject
 * upon request from any calling function. This was designed for 
 * DelayedMatrix objects, to avoid reimplementing arbitrary delayed 
 * R operations in C++; however, it is also useful for unknown matrices.
 */

template<typename T, class V>
class unknown_reader : public dim_checker {
public:    
    unknown_reader(const Rcpp::RObject&);
    ~unknown_reader();

    T get(size_t, size_t);

    template <class Iter> 
    void get_row(size_t, Iter, size_t, size_t);

    template <class Iter>
    void get_col(size_t, Iter, size_t, size_t);

    Rcpp::RObject yield() const;
    matrix_type get_matrix_type () const;
private:
    Rcpp::RObject original;

    Rcpp::Environment beachenv;
    Rcpp::Function realizer_row, realizer_col;
    
    V storage;
    void update_storage_by_row(size_t);
    void update_storage_by_col(size_t);

    Rcpp::IntegerVector row_indices, col_indices;
    int chunk_nrow, chunk_ncol;
};

template<typename T, class V>
unknown_reader<T, V>::unknown_reader(const Rcpp::RObject& in) : original(in), 
    beachenv(Rcpp::Environment::namespace_env("beachmat")),
    realizer_row(beachenv["realizeDelayedMatrixByRow"]), realizer_col(beachenv["realizeDelayedMatrixByCol"]),
    row_indices(2), col_indices(2), chunk_nrow(0), chunk_ncol(0) {

    Rcpp::Function getdims(beachenv["setupDelayedMatrix"]);
    Rcpp::List dimdata=getdims(in);

    Rcpp::IntegerVector matdims(dimdata[0]);
    this->fill_dims(matdims);

    Rcpp::IntegerVector chunkdims(dimdata[1]);
    chunk_nrow=chunkdims[0];
    chunk_ncol=chunkdims[1];
    return;
}

template<typename T, class V>
unknown_reader<T, V>::~unknown_reader() {}

template<typename T, class V>
void unknown_reader<T, V>::update_storage_by_row(size_t r) {
    if (r < row_indices[0] || r >= row_indices[1]) {
        row_indices[0] = std::floor(r/chunk_nrow) * chunk_nrow;
        row_indices[1] = std::min(row_indices[0] + chunk_nrow, int(this->nrow));
        storage=realizer_row(original, row_indices); 
    }
}

template<typename T, class V>
void unknown_reader<T, V>::update_storage_by_col(size_t c) {
    if (c < col_indices[0] || c >= col_indices[1]) {
        col_indices[0] = std::floor(c/chunk_ncol) * chunk_ncol;
        col_indices[1] = std::min(col_indices[0] + chunk_ncol, int(this->ncol));
        storage=realizer_col(original, col_indices); 
    }
}

template<typename T, class V>
T unknown_reader<T, V>::get(size_t r, size_t c) {
    check_oneargs(r, c);
    update_storage_by_col(c);
    return storage[r + this->nrow * (c - size_t(col_indices[0]))]; 
}

template<typename T, class V>
template <class Iter>
void unknown_reader<T, V>::get_row(size_t r, Iter out, size_t first, size_t last) {
    check_rowargs(r, first, last);
    update_storage_by_row(r);
    auto src=storage.begin() + r - size_t(row_indices[0]);
    for (size_t col=first; col<last; ++col, src+=chunk_nrow, ++out) { (*out)=(*src); }
    return;
}
 
template<typename T, class V>
template <class Iter>
void unknown_reader<T, V>::get_col(size_t c, Iter out, size_t first, size_t last) {
    check_colargs(c, first, last);
    update_storage_by_col(c);
    auto src=storage.begin() + (c - size_t(col_indices[0])) * (this->nrow);
    std::copy(src + first, src + last, out);
    return;
}
    
template<typename T, class V>
Rcpp::RObject unknown_reader<T, V>::yield() const {
    return original;
}

template<typename t, class v>
matrix_type unknown_reader<t, v>::get_matrix_type () const {
    return UNKNOWN;
}

}

#endif
