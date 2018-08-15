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

    void get_row(size_t, int*, size_t, size_t);
    void get_row(size_t, double*, size_t, size_t);

    void get_col(size_t, int*, size_t, size_t);
    void get_col(size_t, double*, size_t, size_t);

    void get_rows(int*, size_t, int*, size_t, size_t);
    void get_rows(int*, size_t, double*, size_t, size_t);

    void get_cols(int*, size_t, int*, size_t, size_t);
    void get_cols(int*, size_t, double*, size_t, size_t);

    Rcpp::RObject yield() const;
    matrix_type get_matrix_type () const;
private:
    Rcpp::RObject original;
    void * ptr;

    T (*load) (void *, int, int);
    
    void (*load_col_int) (void *, int, int*, int, int);
    void (*load_row_int) (void *, int, int*, int, int);
    void (*load_col_dbl) (void *, int, double*, int, int);
    void (*load_row_dbl) (void *, int, double*, int, int);

    void (*load_cols_int) (void *, int*, int, int*, int, int);
    void (*load_rows_int) (void *, int*, int, int*, int, int);
    void (*load_cols_dbl) (void *, int*, int, double*, int, int);
    void (*load_rows_dbl) (void *, int*, int, double*, int, int);

    void * (*clone) (void *);
    void (*destroy) (void *);
};

/* Rule of 5 */

template<typename T, class V>
external_lin_reader<T, V>::external_lin_reader(const Rcpp::RObject& incoming) : original(incoming) {
    std::string data_type=translate_type(V(0).sexp_type());
    const char* type=data_type.c_str();

    // Getting the package of origin.
    auto classinfo=get_class_package(incoming);
    const char* pkg=classinfo.second.c_str();

    // Getting all required functions from the corresponding shared library.
    load=reinterpret_cast<T (*)(void *, int, int)>(R_GetCCallable(pkg, combine_strings("load_", type).c_str()));

    load_col_int=reinterpret_cast<void (*)(void *, int, int*, int, int)>(R_GetCCallable(pkg, combine_strings("load_col2int_", type).c_str()));
    load_row_int=reinterpret_cast<void (*)(void *, int, int*, int, int)>(R_GetCCallable(pkg, combine_strings("load_row2int_", type).c_str()));
    load_col_dbl=reinterpret_cast<void (*)(void *, int, double*, int, int)>(R_GetCCallable(pkg, combine_strings("load_col2dbl_", type).c_str()));
    load_row_dbl=reinterpret_cast<void (*)(void *, int, double*, int, int)>(R_GetCCallable(pkg, combine_strings("load_row2dbl_", type).c_str()));

    load_cols_int=reinterpret_cast<void (*)(void *, int*, int, int*, int, int)>(R_GetCCallable(pkg, combine_strings("load_cols2int_", type).c_str()));
    load_rows_int=reinterpret_cast<void (*)(void *, int*, int, int*, int, int)>(R_GetCCallable(pkg, combine_strings("load_rows2int_", type).c_str()));
    load_cols_dbl=reinterpret_cast<void (*)(void *, int*, int, double*, int, int)>(R_GetCCallable(pkg, combine_strings("load_cols2dbl_", type).c_str()));
    load_rows_dbl=reinterpret_cast<void (*)(void *, int*, int, double*, int, int)>(R_GetCCallable(pkg, combine_strings("load_rows2dbl_", type).c_str()));

    clone=reinterpret_cast<void * (*)(void *)>(R_GetCCallable(pkg, combine_strings("clone_", type).c_str()));
    destroy=reinterpret_cast<void (*)(void *)>(R_GetCCallable(pkg, combine_strings("destroy_", type).c_str()));
   
    // Allocating memory last, so we don't have to handle memory deallocation if the above steps throw.
    auto create=reinterpret_cast<void * (*)(SEXP)>(R_GetCCallable(pkg, combine_strings("clone_", type).c_str()));
    ptr=create(original);
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
    destroy{other.destroy}
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
    return load(ptr, r, c);
}

template<typename T, class V>
void external_lin_reader<T, V>::get_row(size_t r, int* out, size_t first, size_t last) {
    load_row_int(ptr, r, out, first, last);
    return;
}

template<typename T, class V>
void external_lin_reader<T, V>::get_row(size_t r, double* out, size_t first, size_t last) {
    load_row_dbl(ptr, r, out, first, last);
    return;
}

template<typename T, class V>
void external_lin_reader<T, V>::get_col(size_t c, int* out, size_t first, size_t last) {
    load_col_int(ptr, c, out, first, last);
    return;
}

template<typename T, class V>
void external_lin_reader<T, V>::get_col(size_t c, double* out, size_t first, size_t last) {
    load_col_dbl(ptr, c, out, first, last);
    return;
}

/* Multi getters. */

template<typename T, class V>
void external_lin_reader<T, V>::get_rows(int* r, size_t n, int* out, size_t first, size_t last) {
    load_rows_int(ptr, r, n, out, first, last);
    return;
}

template<typename T, class V>
void external_lin_reader<T, V>::get_rows(int* r, size_t n, double* out, size_t first, size_t last) {
    load_rows_dbl(ptr, r, n, out, first, last);
    return;
}

template<typename T, class V>
void external_lin_reader<T, V>::get_cols(int* c, size_t n, int* out, size_t first, size_t last) {
    load_cols_int(ptr, c, n, out, first, last);
    return;
}

template<typename T, class V>
void external_lin_reader<T, V>::get_cols(int* c, size_t n, double* out, size_t first, size_t last) {
    load_cols_dbl(ptr, c, n, out, first, last);
    return;
}

/* Miscellaneous functions. */

template<typename T, class V>
Rcpp::RObject external_lin_reader<T, V>::yield() const {
    return original;
}

template<typename T, class V>
matrix_type external_lin_reader<T, V>::get_matrix_type () const {
    return UNKNOWN;
}

}

#endif
