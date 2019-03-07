#ifndef BEACHMAT_EXTERNAL_READER_H
#define BEACHMAT_EXTERNAL_READER_H

#include "Rcpp.h"

#include "../utils/dim_checker.h"
#include "../utils/external.h"

#include <string>

namespace beachmat {

/******************************
 *** Basic class definition ***
 ******************************/

template<typename T, class V>
class external_reader_base : public dim_checker {
public:
    external_reader_base(const Rcpp::RObject& incoming) : original(incoming) {
        auto data_type=get_type();
        const char* type=data_type.c_str();
        auto classinfo=get_class_package(original);
        const char* matclass=classinfo.first.c_str();
        const char* pkg=classinfo.second.c_str();

        // Getting required functions from the corresponding shared library.
        auto load_name=get_external_name(matclass, type, "input", "get");
        load=reinterpret_cast<void (*)(void *, size_t, size_t, T*)>(R_GetCCallable(pkg, load_name.c_str()));

        auto load_const_name=get_external_name(matclass, type, "input", "getConstCol");
        load_const_col=reinterpret_cast<void (*)(void *, size_t, typename V::iterator*, size_t, size_t, typename V::iterator*)>(
            R_GetCCallable(pkg, load_const_name.c_str()));

        auto load_const_indexed_name=get_external_name(matclass, type, "input", "getConstColIndexed");
        load_const_col_indexed=reinterpret_cast<size_t (*)(void *, size_t, Rcpp::IntegerVector::iterator*, typename V::iterator*, size_t, size_t)>(
            R_GetCCallable(pkg, load_const_indexed_name.c_str()));

        ex=external_ptr(original, pkg, matclass, type); // move assignment.

        // Getting the dimensions from the created object.
        auto get_dim_name=get_external_name(matclass, type, "input", "dim");
        auto dimgetter=reinterpret_cast<void (*)(void*, size_t*, size_t*)>(R_GetCCallable(pkg, get_dim_name.c_str()));
        dimgetter(ex.get(), &nrow, &ncol);
        return;
    }

    ~external_reader_base() = default;
    external_reader_base(const external_reader_base&) = default;
    external_reader_base& operator=(const external_reader_base&) = default;
    external_reader_base(external_reader_base&&) = default;
    external_reader_base& operator=(external_reader_base&&) = default;

    T get(size_t r, size_t c) {
        this->check_oneargs(r, c);
        T output;
        load(ex.get(), r, c, &output);
        return output;
    }

    typename V::iterator get_const_col(size_t c, typename V::iterator out, size_t first, size_t last) {
        this->check_colargs(c, first, last);
        typename V::iterator output;
        load_const_col(ex.get(), c, &out, first, last, &output);
        return output;
    }

    size_t get_const_col_indexed(size_t c, Rcpp::IntegerVector::iterator& iIt, typename V::iterator& vIt, size_t first, size_t last) {
        this->check_colargs(c, first, last);
        return load_const_col_indexed(ex.get(), c, &iIt, &vIt, first, last);
    }

    Rcpp::RObject yield() const {
        return original;
    }

    matrix_type get_matrix_type () const {
        return EXTERNAL;
    }

protected:
    Rcpp::RObject original;

    void (*load) (void *, size_t, size_t, T*);
    void (*load_const_col) (void *, size_t, typename V::iterator *, size_t, size_t, typename V::iterator*);
    size_t (*load_const_col_indexed) (void *, size_t, Rcpp::IntegerVector::iterator*, typename V::iterator*, size_t, size_t);

    // Getting the type.
    static std::string get_type() {
        return translate_type(V(0).sexp_type());
    }

protected:
    external_ptr ex;
};

/*******************************
 *** Class with solo getters ***
 *******************************/

template<typename T, class V>
class external_reader : public external_reader_base<T, V> {
private:
    // Typedef'd for convenience.
    typedef Rcpp::IntegerVector::iterator RcppIntIt;
    typedef typename V::iterator RcppValIt;

    void (*load_col) (void *, size_t, RcppValIt*, size_t, size_t);
    void (*load_row) (void *, size_t, RcppValIt*, size_t, size_t);

    void (*load_cols) (void *, RcppIntIt*, size_t, RcppValIt*, size_t, size_t);
    void (*load_rows) (void *, RcppIntIt*, size_t, RcppValIt*, size_t, size_t);

public:    
    external_reader(const Rcpp::RObject& incoming) : external_reader_base<T, V>(incoming) {
        auto data_type=this->get_type();
        const char* type=data_type.c_str();
        auto classinfo=get_class_package(this->original);
        const char* matclass=classinfo.first.c_str();
        const char* pkg=classinfo.second.c_str();
    
        // Getting all required functions from the corresponding shared library.
        auto load_col_name=get_external_name(matclass, type, "input", "getCol");
        load_col=reinterpret_cast<void (*)(void *, size_t, RcppValIt*, size_t, size_t)>(R_GetCCallable(pkg, load_col_name.c_str()));
    
        auto load_row_name=get_external_name(matclass, type, "input", "getRow");
        load_row=reinterpret_cast<void (*)(void *, size_t, RcppValIt*, size_t, size_t)>(R_GetCCallable(pkg, load_row_name.c_str()));
    
        auto load_cols_name=get_external_name(matclass, type, "input", "getCols");
        load_cols=reinterpret_cast<void (*)(void *, RcppIntIt*, size_t, RcppValIt*, size_t, size_t)>(R_GetCCallable(pkg, load_cols_name.c_str()));
    
        auto load_rows_name=get_external_name(matclass, type, "input", "getRows");
        load_rows=reinterpret_cast<void (*)(void *, RcppIntIt*, size_t, RcppValIt*, size_t, size_t)>(R_GetCCallable(pkg, load_rows_name.c_str()));

        return;
    }

    ~external_reader() = default;
    external_reader(const external_reader&) = default;
    external_reader& operator=(const external_reader&) = default;
    external_reader(external_reader&&) = default;
    external_reader& operator=(external_reader&&) = default;

    void get_row(size_t r, RcppValIt out, size_t first, size_t last) {
        this->check_rowargs(r, first, last);
        load_row(this->ex.get(), r, &out, first, last);
        return;
    }

    void get_col(size_t c, RcppValIt out, size_t first, size_t last) {
        this->check_colargs(c, first, last);
        load_col(this->ex.get(), c, &out, first, last);
        return;
    }

    void get_rows(RcppIntIt rIt, size_t n, RcppValIt out, size_t first, size_t last) {
        this->check_rowargs(0, first, last);
        this->check_row_indices(rIt, n);
        load_rows(this->ex.get(), &rIt, n, &out, first, last);
        return;
    }

    void get_cols(RcppIntIt cIt, size_t n, RcppValIt out, size_t first, size_t last) {
        this->check_colargs(0, first, last);
        this->check_col_indices(cIt, n);
        load_cols(this->ex.get(), &cIt, n, &out, first, last);
        return;
    }
};

/******************************
 *** Class with LIN getters ***
 ******************************/

template<typename T, class V>
class external_lin_reader : public external_reader_base<T, V> {
private:
    // Typedef'd for convenience.
    typedef Rcpp::IntegerVector::iterator RcppIntIt;
    typedef Rcpp::NumericVector::iterator RcppNumIt;

    void (*load_col_int) (void *, size_t, RcppIntIt*, size_t, size_t);
    void (*load_row_int) (void *, size_t, RcppIntIt*, size_t, size_t);
    void (*load_col_dbl) (void *, size_t, RcppNumIt*, size_t, size_t);
    void (*load_row_dbl) (void *, size_t, RcppNumIt*, size_t, size_t);

    void (*load_cols_int) (void *, RcppIntIt*, size_t, RcppIntIt*, size_t, size_t);
    void (*load_rows_int) (void *, RcppIntIt*, size_t, RcppIntIt*, size_t, size_t);
    void (*load_cols_dbl) (void *, RcppIntIt*, size_t, RcppNumIt*, size_t, size_t);
    void (*load_rows_dbl) (void *, RcppIntIt*, size_t, RcppNumIt*, size_t, size_t);
public:    
    external_lin_reader(const Rcpp::RObject& incoming) : external_reader_base<T, V>(incoming) {
        auto data_type=this->get_type();
        const char* type=data_type.c_str();
        auto classinfo=get_class_package(this->original);
        const char* matclass=classinfo.first.c_str();
        const char* pkg=classinfo.second.c_str();

        // Getting all required functions from the corresponding shared library.
        auto load_col2int_name=get_external_name(matclass, type, "input", "getCol", "integer");
        load_col_int=reinterpret_cast<void (*)(void *, size_t, RcppIntIt*, size_t, size_t)>(R_GetCCallable(pkg, load_col2int_name.c_str()));

        auto load_row2int_name=get_external_name(matclass, type, "input", "getRow", "integer");
        load_row_int=reinterpret_cast<void (*)(void *, size_t, RcppIntIt*, size_t, size_t)>(R_GetCCallable(pkg, load_row2int_name.c_str()));

        auto load_col2dbl_name=get_external_name(matclass, type, "input", "getCol", "numeric");
        load_col_dbl=reinterpret_cast<void (*)(void *, size_t, RcppNumIt*, size_t, size_t)>(R_GetCCallable(pkg, load_col2dbl_name.c_str()));

        auto load_row2dbl_name=get_external_name(matclass, type, "input", "getRow", "numeric");
        load_row_dbl=reinterpret_cast<void (*)(void *, size_t, RcppNumIt*, size_t, size_t)>(R_GetCCallable(pkg, load_row2dbl_name.c_str()));

        auto load_cols2int_name=get_external_name(matclass, type, "input", "getCols", "integer");
        load_cols_int=reinterpret_cast<void (*)(void *, RcppIntIt*, size_t, RcppIntIt*, size_t, size_t)>(R_GetCCallable(pkg, load_cols2int_name.c_str()));

        auto load_rows2int_name=get_external_name(matclass, type, "input", "getRows", "integer");
        load_rows_int=reinterpret_cast<void (*)(void *, RcppIntIt*, size_t, RcppIntIt*, size_t, size_t)>(R_GetCCallable(pkg, load_rows2int_name.c_str()));

        auto load_cols2dbl_name=get_external_name(matclass, type, "input", "getCols", "numeric");
        load_cols_dbl=reinterpret_cast<void (*)(void *, RcppIntIt*, size_t, RcppNumIt*, size_t, size_t)>(R_GetCCallable(pkg, load_cols2dbl_name.c_str()));

        auto load_rows2dbl_name=get_external_name(matclass, type, "input", "getRows", "numeric");
        load_rows_dbl=reinterpret_cast<void (*)(void *, RcppIntIt*, size_t, RcppNumIt*, size_t, size_t)>(R_GetCCallable(pkg, load_rows2dbl_name.c_str()));

        return;
    }

    ~external_lin_reader() = default;
    external_lin_reader(const external_lin_reader&) = default;
    external_lin_reader& operator=(const external_lin_reader&) = default;
    external_lin_reader(external_lin_reader&&) = default;
    external_lin_reader& operator=(external_lin_reader&&) = default;

    // Basic getters
    void get_row(size_t r, RcppIntIt out, size_t first, size_t last) {
        this->check_rowargs(r, first, last);
        load_row_int(this->ex.get(), r, &out, first, last);
        return;
    }

    void get_row(size_t r, RcppNumIt out, size_t first, size_t last) {
        this->check_rowargs(r, first, last);
        load_row_dbl(this->ex.get(), r, &out, first, last);
        return;
    }

    void get_col(size_t c, RcppIntIt out, size_t first, size_t last) {
        this->check_colargs(c, first, last);
        load_col_int(this->ex.get(), c, &out, first, last);
        return;
    }

    void get_col(size_t c, RcppNumIt out, size_t first, size_t last) {
        this->check_colargs(c, first, last);
        load_col_dbl(this->ex.get(), c, &out, first, last);
        return;
    }

    // Multi getters
    void get_rows(RcppIntIt rIt, size_t n, RcppIntIt out, size_t first, size_t last) {
        this->check_rowargs(0, first, last);
        this->check_row_indices(rIt, n);
        load_rows_int(this->ex.get(), &rIt, n, &out, first, last);
        return;
    }
    
    void get_rows(RcppIntIt rIt, size_t n, RcppNumIt out, size_t first, size_t last) {
        this->check_rowargs(0, first, last);
        this->check_row_indices(rIt, n);
        load_rows_dbl(this->ex.get(), &rIt, n, &out, first, last);
        return;
    }
    
    void get_cols(RcppIntIt cIt, size_t n, RcppIntIt out, size_t first, size_t last) {
        this->check_colargs(0, first, last);
        this->check_col_indices(cIt, n);
        load_cols_int(this->ex.get(), &cIt, n, &out, first, last);
        return;
    }
    
    void get_cols(RcppIntIt cIt, size_t n, RcppNumIt out, size_t first, size_t last) {
        this->check_colargs(0, first, last);
        this->check_col_indices(cIt, n);
        load_cols_dbl(this->ex.get(), &cIt, n, &out, first, last);
        return;
    }
};

}

#endif
