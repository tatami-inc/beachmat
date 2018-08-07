#ifndef UNKNOWN_READER_H
#define UNKNOWN_READER_H

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

    template <class Iter>
    void get_rows(Rcpp::IntegerVector::iterator, size_t, Iter, size_t, size_t);

    template <class Iter>
    void get_cols(Rcpp::IntegerVector::iterator, size_t, Iter, size_t, size_t);

    Rcpp::RObject yield() const;
    matrix_type get_matrix_type () const;
private:
    Rcpp::RObject original;
    Rcpp::Environment beachenv;
    Rcpp::Function realizer;
    
    V storage; // this does not need to be copyable, as the entire object is replaced when a new block is loaded.
    void update_storage_by_row(size_t, size_t, size_t);
    void update_storage_by_col(size_t, size_t, size_t);
    size_t get_row_storage_start_row() const;
    size_t get_row_storage_start_col() const;
    size_t get_row_storage_ncols() const;
    size_t get_col_storage_start_row() const;
    size_t get_col_storage_start_col() const;
    size_t get_col_storage_nrows() const;

    copyable_holder<Rcpp::IntegerVector> row_indices, row_slices;
    copyable_holder<Rcpp::IntegerVector> col_indices, col_slices;
    copyable_holder<Rcpp::LogicalVector> do_transpose;
    int chunk_nrow, chunk_ncol;
};

template<typename T, class V>
unknown_reader<T, V>::unknown_reader(const Rcpp::RObject& in) : original(in), 
        beachenv(Rcpp::Environment::namespace_env("beachmat")), realizer(beachenv["realizeByRange"]), 
        row_indices(2), row_slices(2), col_indices(2), col_slices(2), do_transpose(1), 
        chunk_nrow(0), chunk_ncol(0) {

    Rcpp::Function getdims(beachenv["setupUnknownMatrix"]);
    Rcpp::List dimdata=getdims(in);

    Rcpp::IntegerVector matdims(dimdata[0]);
    this->fill_dims(matdims);

    Rcpp::IntegerVector chunkdims(dimdata[1]);
    chunk_nrow=chunkdims[0];
    chunk_ncol=chunkdims[1];

    do_transpose.vec[0]=1;
    return;
}

template<typename T, class V>
unknown_reader<T, V>::~unknown_reader() {}

/* Define storage-related methods. */

template<typename T, class V>
void unknown_reader<T, V>::update_storage_by_row(size_t r, size_t first, size_t last) {
    auto& rinds=row_indices.vec;
    auto& cinds=col_slices.vec;

    if (r < size_t(rinds[0]) || r >= size_t(rinds[1]) || first < size_t(cinds[0]) || last > size_t(cinds[1])) {
        rinds[0] = std::floor(r/chunk_nrow) * chunk_nrow;
        rinds[1] = std::min(chunk_nrow, int(this->nrow) - rinds[0]);

        cinds[0] = first;
        cinds[1] = last - first;
        storage = realizer(original, rinds, cinds, do_transpose.vec); // Transposed, so storage is effectively row-major!

        // Resetting for easier checks above.
        rinds[1] += rinds[0];
        cinds[1] = last;
    }
    return;
}

template<typename T, class V>
void unknown_reader<T, V>::update_storage_by_col(size_t c, size_t first, size_t last) {
    auto& cinds=col_indices.vec;
    auto& rinds=row_slices.vec;

    if (c < size_t(cinds[0]) || c >= size_t(cinds[1]) || first < size_t(rinds[0]) || last > size_t(rinds[1])) {
        cinds[0] = std::floor(c/chunk_ncol) * chunk_ncol;
        cinds[1] = std::min(chunk_ncol, int(this->ncol) - cinds[0]);

        rinds[0] = first;
        rinds[1] = last - first;
        storage = realizer(original, rinds, cinds);

        // Resetting for easier checks above.
        cinds[1] += cinds[0];
        rinds[1] = last;
    }
    return;
}

template<typename T, class V>
size_t unknown_reader<T, V>::get_row_storage_start_row() const { return row_indices.vec[0]; }

template<typename T, class V>
size_t unknown_reader<T, V>::get_row_storage_start_col() const { return col_slices.vec[0]; }

template<typename T, class V>
size_t unknown_reader<T, V>::get_row_storage_ncols() const { return col_slices.vec[1] - col_slices.vec[0]; }

template<typename T, class V>
size_t unknown_reader<T, V>::get_col_storage_start_row() const { return row_slices.vec[0]; }

template<typename T, class V>
size_t unknown_reader<T, V>::get_col_storage_start_col() const { return col_indices.vec[0]; }

template<typename T, class V>
size_t unknown_reader<T, V>::get_col_storage_nrows() const { return row_slices.vec[1] - row_slices.vec[0]; }

/*** Basic getter methods ***/

template<typename T, class V>
T unknown_reader<T, V>::get(size_t r, size_t c) {
    check_oneargs(r, c);
    update_storage_by_col(c, 0, this->nrow); // keeps the whole block for further 'get()' queries.
    return storage[(c - this->get_col_storage_start_col()) * this->nrow + r];
}

template<typename T, class V>
template <class Iter>
void unknown_reader<T, V>::get_row(size_t r, Iter out, size_t first, size_t last) {
    check_rowargs(r, first, last);
    update_storage_by_row(r, first, last);

    // It's effectively row-major storage due the transposition.
    auto src=storage.begin() + 
        (r - this->get_row_storage_start_row()) * (this->get_row_storage_ncols()) +
        (first - this->get_row_storage_start_col());
    std::copy(src, src + (last - first), out);
    return;
}
 
template<typename T, class V>
template <class Iter>
void unknown_reader<T, V>::get_col(size_t c, Iter out, size_t first, size_t last) {
    check_colargs(c, first, last);
    update_storage_by_col(c, first, last);

    auto src=storage.begin() + 
        (c - this->get_col_storage_start_col()) * (this->get_col_storage_nrows()) +
        (first - this->get_col_storage_start_row());
    std::copy(src, src + (last - first), out);
    return;
}

/*** Multi getter methods ***/

template<typename T, class V>
template<class Iter>
void unknown_reader<T, V>::get_rows(Rcpp::IntegerVector::iterator cIt, size_t n, Iter out, size_t first, size_t last) {
    check_rowargs(0, first, last);
    check_row_indices(cIt, n);

    // Need to make a copy of the indexed (1-indexed) to pass to the function.
    // Don't use col_indices.vec, avoid bugs upon interaction with update_storage().
    Rcpp::IntegerVector cur_indices(cIt, cIt+n);
    for (auto& i : cur_indices) { ++i; }
    Rcpp::IntegerVector col_range=Rcpp::IntegerVector::create(first, last-first);
    
    Rcpp::Function indexed_realizer(beachenv["realizeByIndexRange"]);
    V tmp_store=indexed_realizer(original, cur_indices, col_range);
    std::copy(tmp_store.begin(), tmp_store.end(), out);
    return;
}

template<typename T, class V>
template<class Iter>
void unknown_reader<T, V>::get_cols(Rcpp::IntegerVector::iterator cIt, size_t n, Iter out, size_t first, size_t last) {
    check_colargs(0, first, last);
    check_col_indices(cIt, n);

    // Need to make a copy of the indices (1-indexed) to pass to the function.
    // Don't use row_indices.vec, avoid bugs upon interaction with update_storage().
    Rcpp::IntegerVector cur_indices(cIt, cIt+n);
    for (auto& i : cur_indices) { ++i; }
    Rcpp::IntegerVector row_range=Rcpp::IntegerVector::create(first, last-first);
    
    Rcpp::Function indexed_realizer(beachenv["realizeByRangeIndex"]);
    V tmp_store=indexed_realizer(original, row_range, cur_indices);
    std::copy(tmp_store.begin(), tmp_store.end(), out);
    return;
}

/*** Other methods ***/

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
