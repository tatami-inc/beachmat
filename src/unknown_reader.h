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
    ~unknown_reader() = default;
    unknown_reader(const unknown_reader&) = default;
    unknown_reader& operator=(const unknown_reader&) = default;
    unknown_reader(unknown_reader&&) = default;
    unknown_reader& operator=(unknown_reader&&) = default;

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
    
    size_t storage_start_row, storage_end_row, storage_start_col, storage_end_col; 
    bool oncol;

    static size_t get_upper_chunk_bound (size_t end, size_t chunk_dim);

    copyable_holder<Rcpp::IntegerVector> indices, slices;
    copyable_holder<Rcpp::LogicalVector> do_transpose;
    size_t chunk_nrow, chunk_ncol;
};

/* Constructor definition. */

template<typename T, class V>
unknown_reader<T, V>::unknown_reader(const Rcpp::RObject& in) : original(in), 
        beachenv(Rcpp::Environment::namespace_env("beachmat")), realizer(beachenv["realizeByRange"]), 
        storage_start_row(0), storage_end_row(0), storage_start_col(0), storage_end_col(0), 
        oncol(false), indices(2), slices(2), do_transpose(1), 
        chunk_nrow(0), chunk_ncol(0) {

    Rcpp::Function getdims(beachenv["setupUnknownMatrix"]);
    Rcpp::List dimdata=getdims(in);

    Rcpp::IntegerVector matdims(dimdata[0]);
    this->fill_dims(matdims);

    Rcpp::IntegerVector chunkdims(dimdata[1]);
    chunk_nrow=std::max(1, chunkdims[0]);
    chunk_ncol=std::max(1, chunkdims[1]);

    do_transpose.vec[0]=1;
    return;
}

/* Define storage-related methods. */

template<typename T, class V>
size_t unknown_reader<T, V>::get_upper_chunk_bound (size_t end, size_t chunk_dim) {
    if (end <= chunk_dim) {
        return chunk_dim;
    }
    return (static_cast<size_t>((end - 1)/chunk_dim) + 1) * chunk_dim; // implicit and intended floor upon division.
}

template<typename T, class V>
void unknown_reader<T, V>::update_storage_by_row(size_t r, size_t first, size_t last) {
    // We assume that all inputs are valid, as check_* functions are called elsewhere.
    if (oncol 
            || r < storage_start_row || r >= storage_end_row 
            || first < storage_start_col || last > storage_end_col) { 

        storage_start_row = static_cast<size_t>(r/chunk_nrow) * chunk_nrow; // implicit floor in the division.
        storage_end_row = std::min(chunk_nrow + storage_start_row, this->nrow);

        indices.vec[0] = storage_start_row; 
        indices.vec[1] = storage_end_row - storage_start_row;

        storage_start_col = static_cast<size_t>(first/chunk_ncol) * chunk_ncol;
        storage_end_col = std::min(get_upper_chunk_bound(last, chunk_ncol), this->ncol);

        slices.vec[0] = storage_start_col;
        slices.vec[1] = storage_end_col - storage_start_col;

        storage = realizer(original, indices.vec, slices.vec, do_transpose.vec); // Transposed, so storage is effectively row-major!
        oncol=false;
    }
    return;
}

template<typename T, class V>
void unknown_reader<T, V>::update_storage_by_col(size_t c, size_t first, size_t last) {
    if (!oncol 
            || c < storage_start_col || c >= storage_end_col 
            || first < storage_start_row || last > storage_end_row) { 
    
        storage_start_col = static_cast<size_t>(c/chunk_ncol) * chunk_ncol; // implicit floor in the division. 
        storage_end_col = std::min(chunk_ncol + storage_start_col, this->ncol);

        indices.vec[0] = storage_start_col;
        indices.vec[1] = storage_end_col - storage_start_col;

        storage_start_row = static_cast<size_t>(first/chunk_nrow) * chunk_nrow;
        storage_end_row = std::min(get_upper_chunk_bound(last, chunk_nrow), this->nrow);

        slices.vec[0] = storage_start_row; 
        slices.vec[1] = storage_end_row - storage_start_row;

        storage = realizer(original, slices.vec, indices.vec);
        oncol=true;
    }
    return;
}

/*** Basic getter methods ***/

template<typename T, class V>
T unknown_reader<T, V>::get(size_t r, size_t c) {
    check_oneargs(r, c);
    update_storage_by_col(c, 0, this->nrow); // keeps the whole block for further 'get()' queries.
    return storage[(c - storage_start_col) * this->nrow + r];
}

template<typename T, class V>
template <class Iter>
void unknown_reader<T, V>::get_row(size_t r, Iter out, size_t first, size_t last) {
    check_rowargs(r, first, last);
    update_storage_by_row(r, first, last);

    // It's effectively row-major storage due the transposition.
    auto src=storage.begin() + 
        (r - storage_start_row) * (storage_end_col - storage_start_col) +
        (first - storage_start_col);
    std::copy(src, src + (last - first), out);
    return;
}
 
template<typename T, class V>
template <class Iter>
void unknown_reader<T, V>::get_col(size_t c, Iter out, size_t first, size_t last) {
    check_colargs(c, first, last);
    update_storage_by_col(c, first, last);

    auto src=storage.begin() + 
        (c - storage_start_col) * (storage_end_row - storage_start_row) + 
        (first - storage_start_row);
    std::copy(src, src + (last - first), out);
    return;
}

/*** Multi getter methods ***/

template<typename T, class V>
template<class Iter>
void unknown_reader<T, V>::get_rows(Rcpp::IntegerVector::iterator rIt, size_t n, Iter out, size_t first, size_t last) {
    check_rowargs(0, first, last);
    check_row_indices(rIt, n);

    // Need to make a copy of the indexed (1-indexed) to pass to the function.
    Rcpp::IntegerVector cur_indices(rIt, rIt+n);
    for (auto& i : cur_indices) { ++i; }
    slices.vec[0]=first;
    slices.vec[1]=last-first;
    
    Rcpp::Function indexed_realizer(beachenv["realizeByIndexRange"]);
    V tmp_store=indexed_realizer(original, cur_indices, slices.vec);
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
    slices.vec[0]=first;
    slices.vec[1]=last-first;

    Rcpp::Function indexed_realizer(beachenv["realizeByRangeIndex"]);
    V tmp_store=indexed_realizer(original, slices.vec, cur_indices);
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
