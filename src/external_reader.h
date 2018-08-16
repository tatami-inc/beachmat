#ifndef EXTERNAL_READER_H
#define EXTERNAL_READER_H

namespace beachmat {

/*** Class definition ***/

template<typename T, class V>
class external_lin_reader : public dim_checker {
public:    
    external_lin_reader(const Rcpp::RObject&);
    ~external_lin_reader();
    external_lin_reader(const external_lin_reader&);
    external_lin_reader& operator=(const external_lin_reader&);
    external_lin_reader(external_lin_reader&&) = default;             // Move constructor
    external_lin_reader& operator=(external_lin_reader&&) = default;  // Move assignment constructor

    T get(size_t, size_t);

    void get_row(size_t, Rcpp::IntegerVector::iterator, size_t, size_t);
    void get_row(size_t, Rcpp::NumericVector::iterator, size_t, size_t);

    void get_col(size_t, Rcpp::IntegerVector::iterator, size_t, size_t);
    void get_col(size_t, Rcpp::NumericVector::iterator, size_t, size_t);

    void get_rows(Rcpp::IntegerVector::iterator, size_t, Rcpp::IntegerVector::iterator, size_t, size_t);
    void get_rows(Rcpp::IntegerVector::iterator, size_t, Rcpp::NumericVector::iterator, size_t, size_t);

    void get_cols(Rcpp::IntegerVector::iterator, size_t, Rcpp::IntegerVector::iterator, size_t, size_t);
    void get_cols(Rcpp::IntegerVector::iterator, size_t, Rcpp::NumericVector::iterator, size_t, size_t);

    Rcpp::RObject yield() const;
    matrix_type get_matrix_type () const;
private:
    Rcpp::RObject original;
    void * ptr;

    T (*load) (void *, size_t, size_t);

    // Typedef'd for convenience.
    typedef Rcpp::IntegerVector::iterator RcppIntIt;
    typedef Rcpp::NumericVector::iterator RcppNumIt;

    void (*load_col_int) (void *, size_t, RcppIntIt, size_t, size_t);
    void (*load_row_int) (void *, size_t, RcppIntIt, size_t, size_t);
    void (*load_col_dbl) (void *, size_t, RcppNumIt, size_t, size_t);
    void (*load_row_dbl) (void *, size_t, RcppNumIt, size_t, size_t);

    void (*load_cols_int) (void *, RcppIntIt, size_t, RcppIntIt, size_t, size_t);
    void (*load_rows_int) (void *, RcppIntIt, size_t, RcppIntIt, size_t, size_t);
    void (*load_cols_dbl) (void *, RcppIntIt, size_t, RcppNumIt, size_t, size_t);
    void (*load_rows_dbl) (void *, RcppIntIt, size_t, RcppNumIt, size_t, size_t);

    void * (*clone) (void *);
    void (*destroy) (void *);
};

/* Rule of 5 */

template<typename T, class V>
external_lin_reader<T, V>::external_lin_reader(const Rcpp::RObject& incoming) : original(incoming) {
    // Getting the current data type.
    std::string data_type=translate_type(V(0).sexp_type());
    const char* type=data_type.c_str();

    // Getting the package of origin.
    auto classinfo=get_class_package(incoming);
    const char* pkg=classinfo.second.c_str();

    // Getting all required functions from the corresponding shared library.
    load=reinterpret_cast<T (*)(void *, size_t, size_t)>(R_GetCCallable(pkg, combine_strings("load_", type).c_str()));

    load_col_int=reinterpret_cast<void (*)(void *, size_t, RcppIntIt, size_t, size_t)>(R_GetCCallable(pkg, combine_strings("load_col2int_", type).c_str()));
    load_row_int=reinterpret_cast<void (*)(void *, size_t, RcppIntIt, size_t, size_t)>(R_GetCCallable(pkg, combine_strings("load_row2int_", type).c_str()));
    load_col_dbl=reinterpret_cast<void (*)(void *, size_t, RcppNumIt, size_t, size_t)>(R_GetCCallable(pkg, combine_strings("load_col2dbl_", type).c_str()));
    load_row_dbl=reinterpret_cast<void (*)(void *, size_t, RcppNumIt, size_t, size_t)>(R_GetCCallable(pkg, combine_strings("load_row2dbl_", type).c_str()));

    load_cols_int=reinterpret_cast<void (*)(void *, RcppIntIt, size_t, RcppIntIt, size_t, size_t)>(R_GetCCallable(pkg, combine_strings("load_cols2int_", type).c_str()));
    load_rows_int=reinterpret_cast<void (*)(void *, RcppIntIt, size_t, RcppIntIt, size_t, size_t)>(R_GetCCallable(pkg, combine_strings("load_rows2int_", type).c_str()));
    load_cols_dbl=reinterpret_cast<void (*)(void *, RcppIntIt, size_t, RcppNumIt, size_t, size_t)>(R_GetCCallable(pkg, combine_strings("load_cols2dbl_", type).c_str()));
    load_rows_dbl=reinterpret_cast<void (*)(void *, RcppIntIt, size_t, RcppNumIt, size_t, size_t)>(R_GetCCallable(pkg, combine_strings("load_rows2dbl_", type).c_str()));

    clone=reinterpret_cast<void * (*)(void *)>(R_GetCCallable(pkg, combine_strings("clone_", type).c_str()));
    destroy=reinterpret_cast<void (*)(void *)>(R_GetCCallable(pkg, combine_strings("destroy_", type).c_str()));
   
    // Allocating memory as late as possible, to minimize the code in the try/catch block.
    auto create=reinterpret_cast<void * (*)(SEXP)>(R_GetCCallable(pkg, combine_strings("create_", type).c_str()));
    ptr=create(original);

    try {
        // Getting the dimensions from the created object.
        auto dimgetter=reinterpret_cast<void (*)(void*, size_t*, size_t*)>(R_GetCCallable(pkg, combine_strings("get_dim_", type).c_str()));
        dimgetter(ptr, &nrow, &ncol);
    } catch (std::exception& e) {
        destroy(ptr);
        throw;
    }

    return;
}

template<typename T, class V>
external_lin_reader<T, V>::~external_lin_reader() {
    destroy(ptr);
    return;
}

template<typename T, class V>
external_lin_reader<T, V>::external_lin_reader(const external_lin_reader& other) : 
    original(other.original),
    ptr(other.clone(other.ptr)),

    load(other.load), 
    load_col_int(other.load_col_int), 
    load_row_int(other.load_row_int), 
    load_col_dbl(other.load_col_dbl), 
    load_row_dbl(other.load_row_dbl),

    load_cols_int(other.load_cols_int), 
    load_rows_int(other.load_rows_int), 
    load_cols_dbl(other.load_cols_dbl), 
    load_rows_dbl(other.load_rows_dbl), 

    clone(other.clone),
    destroy(other.destroy)
{}

template<typename T, class V>
external_lin_reader<T, V>& external_lin_reader<T, V>::operator=(const external_lin_reader& other) { 
    original=other.original;
    ptr=other.clone(other.ptr);

    load=other.load;
    load_col_int=other.load_col_int; 
    load_row_int=other.load_row_int;
    load_col_dbl=other.load_col_dbl; 
    load_row_dbl=other.load_row_dbl;

    load_cols_int=other.load_cols_int; 
    load_rows_int=other.load_rows_int;
    load_cols_dbl=other.load_cols_dbl; 
    load_rows_dbl=other.load_rows_dbl;

    clone=other.clone;
    destroy=other.destroy;
    return *this;
}

/* Basic getters. */

template<typename T, class V>
T external_lin_reader<T, V>::get(size_t r, size_t c) {
    check_oneargs(r, c);
    return load(ptr, r, c);
}

template<typename T, class V>
void external_lin_reader<T, V>::get_row(size_t r, RcppIntIt out, size_t first, size_t last) {
    check_rowargs(r, first, last);
    load_row_int(ptr, r, out, first, last);
    return;
}

template<typename T, class V>
void external_lin_reader<T, V>::get_row(size_t r, RcppNumIt out, size_t first, size_t last) {
    check_rowargs(r, first, last);
    load_row_dbl(ptr, r, out, first, last);
    return;
}

template<typename T, class V>
void external_lin_reader<T, V>::get_col(size_t c, RcppIntIt out, size_t first, size_t last) {
    check_colargs(c, first, last);
    load_col_int(ptr, c, out, first, last);
    return;
}

template<typename T, class V>
void external_lin_reader<T, V>::get_col(size_t c, RcppNumIt out, size_t first, size_t last) {
    check_colargs(c, first, last);
    load_col_dbl(ptr, c, out, first, last);
    return;
}

/* Multi getters. */

template<typename T, class V>
void external_lin_reader<T, V>::get_rows(RcppIntIt rIt, size_t n, RcppIntIt out, size_t first, size_t last) {
    check_rowargs(0, first, last);
    check_row_indices(rIt, n);
    load_rows_int(ptr, rIt, n, out, first, last);
    return;
}

template<typename T, class V>
void external_lin_reader<T, V>::get_rows(RcppIntIt rIt, size_t n, RcppNumIt out, size_t first, size_t last) {
    check_rowargs(0, first, last);
    check_row_indices(rIt, n);
    load_rows_dbl(ptr, rIt, n, out, first, last);
    return;
}

template<typename T, class V>
void external_lin_reader<T, V>::get_cols(RcppIntIt cIt, size_t n, RcppIntIt out, size_t first, size_t last) {
    check_colargs(0, first, last);
    check_col_indices(cIt, n);
    load_cols_int(ptr, cIt, n, out, first, last);
    return;
}

template<typename T, class V>
void external_lin_reader<T, V>::get_cols(RcppIntIt cIt, size_t n, RcppNumIt out, size_t first, size_t last) {
    check_colargs(0, first, last);
    check_col_indices(cIt, n);
    load_cols_dbl(ptr, cIt, n, out, first, last);
    return;
}

/* Miscellaneous functions. */

template<typename T, class V>
Rcpp::RObject external_lin_reader<T, V>::yield() const {
    return original;
}

template<typename T, class V>
matrix_type external_lin_reader<T, V>::get_matrix_type () const {
    return EXTERNAL;
}

/************************
 *** Class definition ***
 ************************/

template<typename T, class V>
class external_reader : public dim_checker {
public:    
    external_reader(const Rcpp::RObject&);
    ~external_reader();
    external_reader(const external_reader&);
    external_reader& operator=(const external_reader&);
    external_reader(external_reader&&) = default;             // Move constructor
    external_reader& operator=(external_reader&&) = default;  // Move assignment constructor

    T get(size_t, size_t);

    void get_row(size_t, Rcpp::StringVector::iterator, size_t, size_t);
    void get_col(size_t, Rcpp::StringVector::iterator, size_t, size_t);

    void get_rows(Rcpp::IntegerVector::iterator, size_t, Rcpp::StringVector::iterator, size_t, size_t);
    void get_cols(Rcpp::IntegerVector::iterator, size_t, Rcpp::StringVector::iterator, size_t, size_t);

    Rcpp::RObject yield() const;
    matrix_type get_matrix_type () const;
private:
    Rcpp::RObject original;
    void * ptr;

    T (*load) (void *, size_t, size_t);

    // Typedef'd for convenience.
    typedef Rcpp::IntegerVector::iterator RcppIntIt;
    typedef Rcpp::StringVector::iterator RcppStrIt;

    void (*load_col_str) (void *, size_t, RcppStrIt, size_t, size_t);
    void (*load_row_str) (void *, size_t, RcppStrIt, size_t, size_t);

    void (*load_cols_str) (void *, RcppIntIt, size_t, RcppStrIt, size_t, size_t);
    void (*load_rows_str) (void *, RcppIntIt, size_t, RcppStrIt, size_t, size_t);

    void * (*clone) (void *);
    void (*destroy) (void *);
};

/* Rule of 5 */

template<typename T, class V>
external_reader<T, V>::external_reader(const Rcpp::RObject& incoming) : original(incoming) {
    // Getting the current data type.
    std::string data_type=translate_type(V(0).sexp_type());
    const char* type=data_type.c_str();

    // Getting the package of origin.
    auto classinfo=get_class_package(incoming);
    const char* pkg=classinfo.second.c_str();

    // Getting all required functions from the corresponding shared library.
    load=reinterpret_cast<T (*)(void *, size_t, size_t)>(R_GetCCallable(pkg, combine_strings("load_", type).c_str()));

    load_col_str=reinterpret_cast<void (*)(void *, size_t, RcppStrIt, size_t, size_t)>(R_GetCCallable(pkg, combine_strings("load_col_", type).c_str()));
    load_row_str=reinterpret_cast<void (*)(void *, size_t, RcppStrIt, size_t, size_t)>(R_GetCCallable(pkg, combine_strings("load_row_", type).c_str()));

    load_cols_str=reinterpret_cast<void (*)(void *, RcppIntIt, size_t, RcppStrIt, size_t, size_t)>(R_GetCCallable(pkg, combine_strings("load_cols_", type).c_str()));
    load_rows_str=reinterpret_cast<void (*)(void *, RcppIntIt, size_t, RcppStrIt, size_t, size_t)>(R_GetCCallable(pkg, combine_strings("load_rows_", type).c_str()));

    clone=reinterpret_cast<void * (*)(void *)>(R_GetCCallable(pkg, combine_strings("clone_", type).c_str()));
    destroy=reinterpret_cast<void (*)(void *)>(R_GetCCallable(pkg, combine_strings("destroy_", type).c_str()));
   
    // Allocating memory as late as possible, to minimize the code in the try/catch block.
    auto create=reinterpret_cast<void * (*)(SEXP)>(R_GetCCallable(pkg, combine_strings("create_", type).c_str()));
    ptr=create(original);

    try {
        // Getting the dimensions from the created object.
        auto dimgetter=reinterpret_cast<void (*)(void*, size_t*, size_t*)>(R_GetCCallable(pkg, combine_strings("get_dim_", type).c_str()));
        dimgetter(ptr, &nrow, &ncol);
    } catch (std::exception& e) {
        destroy(ptr);
        throw;
    }

    return;
}

template<typename T, class V>
external_reader<T, V>::~external_reader() {
    destroy(ptr);
    return;
}

template<typename T, class V>
external_reader<T, V>::external_reader(const external_reader& other) : 
    original(other.original),
    ptr(other.clone(other.ptr)),

    load(other.load), 
    load_col_str(other.load_col_str), 
    load_row_str(other.load_row_str), 

    load_cols_str(other.load_cols_str), 
    load_rows_str(other.load_rows_str), 

    clone(other.clone),
    destroy(other.destroy)
{}

template<typename T, class V>
external_reader<T, V>& external_reader<T, V>::operator=(const external_reader& other) { 
    original=other.original;
    ptr=other.clone(other.ptr);

    load=other.load;
    load_col_str=other.load_col_str; 
    load_row_str=other.load_row_str;

    load_cols_str=other.load_cols_str; 
    load_rows_str=other.load_rows_str;

    clone=other.clone;
    destroy=other.destroy;
    return *this;
}

/* Basic getters. */

template<typename T, class V>
T external_reader<T, V>::get(size_t r, size_t c) {
    check_oneargs(r, c);
    return load(ptr, r, c);
}

template<typename T, class V>
void external_reader<T, V>::get_row(size_t r, RcppStrIt out, size_t first, size_t last) {
    check_rowargs(r, first, last);
    load_row_str(ptr, r, out, first, last);
    return;
}

template<typename T, class V>
void external_reader<T, V>::get_col(size_t c, RcppStrIt out, size_t first, size_t last) {
    check_colargs(c, first, last);
    load_col_str(ptr, c, out, first, last);
    return;
}

/* Multi getters. */

template<typename T, class V>
void external_reader<T, V>::get_rows(RcppIntIt rIt, size_t n, RcppStrIt out, size_t first, size_t last) {
    check_rowargs(0, first, last);
    check_row_indices(rIt, n);
    load_rows_str(ptr, rIt, n, out, first, last);
    return;
}

template<typename T, class V>
void external_reader<T, V>::get_cols(RcppIntIt cIt, size_t n, RcppStrIt out, size_t first, size_t last) {
    check_colargs(0, first, last);
    check_col_indices(cIt, n);
    load_cols_str(ptr, cIt, n, out, first, last);
    return;
}

/* Miscellaneous functions. */

template<typename T, class V>
Rcpp::RObject external_reader<T, V>::yield() const {
    return original;
}

template<typename T, class V>
matrix_type external_reader<T, V>::get_matrix_type () const {
    return EXTERNAL;
}

}

#endif
