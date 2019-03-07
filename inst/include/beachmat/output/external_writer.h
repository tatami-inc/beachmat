#ifndef BEACHMAT_EXTERNAL_WRITER_H
#define BEACHMAT_EXTERNAL_WRITER_H

#include "Rcpp.h"

#include "../utils/utils.h"
#include "../utils/dim_checker.h"
#include "../utils/external.h"

#include <algorithm>

namespace beachmat {

/*************************
 * Base class definition *
 *************************/

template<typename T, class V>
class external_writer_base : public dim_checker {
protected:
    const char * matclass, * type;
    external_ptr ex; 

    void (*store) (void *, size_t, size_t, T*)=NULL;
    void (*load) (void *, size_t, size_t)=NULL;
    void (*report) (void *)=NULL;

public:
    external_writer_base(size_t nr, size_t nc, const char* Pkg, const char* Class, const char* Type) : 
            dim_checker(nr, nc), matclass(Class), type(Type), ex(nr, nc, Pkg, Class, Type) {
        // Define all remaining function pointers.
        auto store_name=get_external_name(matclass, type, "output", "set");
        store=reinterpret_cast<void (*)(void *, size_t, size_t, T*)>(R_GetCCallable(pkg, store_name.c_str()));

        auto load_name=get_external_name(matclass, type, "output", "get");
        load=reinterpret_cast<void (*)(void *, size_t, size_t, T*)>(R_GetCCallable(pkg, load_name.c_str()));

        auto report_name=get_external_name(matclass, type, "output", "yield");
        report_matrix=reinterpret_cast<void (*)(void *)>(R_GetCCallable(pkg, report_name.c_str()));
        return;
    }
    ~external_writer_base() = default;
    external_writer_base(const external_writer_base&) = default;
    external_writer_base& operator=(const external_writer_base&) = default;
    external_writer_base(external_writer_base&&) = default;
    external_writer_base& operator=(external_writer_base&&) = default;

    // Setters:
    void set(size_t r, size_t c, T val) {
        store(ex.get(), r, c, &val);
        return;
    }

    // Getters:
    T get(size_t r, size_t c) {
        return load(ex.get(), r, c);
    }
    
    // Other:
    Rcpp::RObject yield() {
        return yield_matrix(ex.get());
    }

    matrix_type get_matrix_type() const {
        return EXTERNAL;
    }
};

/*******************************
 *** Class with solo getters ***
 *******************************/

template<typename T, class V>
class external_writer : public external_writer_base<T, V> {
public:    
    external_writer(const Rcpp::RObject& incoming) : external_writer_base<T, V>(incoming) {
        auto data_type=this->get_type();
        const char* type=data_type.c_str();
        auto classinfo=get_class_package(this->original);
        const char* pkg=classinfo.second.c_str();
    
        // Getting all required functions from the corresponding shared library.
        auto store_col_name=get_external_name(classinfo, type, "output", "setCol");
        store_col=reinterpret_cast<void (*)(void *, size_t, RcppValIt*, size_t, size_t)>(R_GetCCallable(pkg, store_col_name.c_str()));
    
        auto store_row_name=get_external_name(classinfo, type, "output", "setRow");
        store_row=reinterpret_cast<void (*)(void *, size_t, RcppValIt*, size_t, size_t)>(R_GetCCallable(pkg, store_row_name.c_str()));

        auto store_col_name=get_external_name(classinfo, type, "output", "setColIndexed");
        store_col=reinterpret_cast<void (*)(void *, size_t, size_t, RcppIntIt*, RcppValIt*)>(R_GetCCallable(pkg, store_col_name.c_str()));
    
        auto store_row_name=get_external_name(classinfo, type, "output", "setRowIndexed");
        store_row=reinterpret_cast<void (*)(void *, size_t, size_t, RcppIntIt*, RcppValIt*)>(R_GetCCallable(pkg, store_row_name.c_str()));

        auto load_col_name=get_external_name(classinfo, type, "output", "getCol");
        load_col=reinterpret_cast<void (*)(void *, size_t, RcppValIt*, size_t, size_t)>(R_GetCCallable(pkg, load_col_name.c_str()));
    
        auto load_row_name=get_external_name(classinfo, type, "output", "getRow");
        load_row=reinterpret_cast<void (*)(void *, size_t, RcppValIt*, size_t, size_t)>(R_GetCCallable(pkg, load_row_name.c_str()));

        return;
    }

    ~external_writer() = default;
    external_writer(const external_writer&) = default;
    external_writer& operator=(const external_writer&) = default;
    external_writer(external_writer&&) = default;
    external_writer& operator=(external_writer&&) = default;

    void set_row(size_t r, RcppValIt out, size_t first, size_t last) {
        this->check_rowargs(r, first, last);
        store_row(this->ex.set(), r, &out, first, last);
        return;
    }

    void set_col(size_t c, RcppValIt out, size_t first, size_t last) {
        this->check_colargs(c, first, last);
        store_col(this->ex.set(), c, &out, first, last);
        return;
    }

    void set_col_indexed(size_t c, size_t n, RcppIntIt idx, RcppValIt in) {
        this->check_colargs(c);
        store_col_indexed(ex.get(), c, n, &idx, &in);
        return;
    }

    void set_row_indexed(size_t c, size_t n, RcppIntIt idx, RcppValIt in) {
        this->check_rowargs(c);
        store_row_indexed(ex.get(), c, n, &idx, &in);
        return;
    }

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

private:
    // Typedef'd for convenience.
    typedef Rcpp::IntegerVector::iterator RcppIntIt;
    typedef typename V::iterator RcppValIt;

    void (*store_col) (void *, size_t, RcppValIt*, size_t, size_t);
    void (*store_row) (void *, size_t, RcppValIt*, size_t, size_t);

    void (*store_col_indexed) (void *, size_t, size_t, RcppIntIt*, RcppValIt*);
    void (*store_row_indexed) (void *, size_t, size_t, RcppIntIt*, RcppValIt*);

    void (*load_col) (void *, size_t, RcppValIt*, size_t, size_t);
    void (*load_row) (void *, size_t, RcppValIt*, size_t, size_t);
};

/******************************
 *** Class with LIN getters ***
 ******************************/

template<typename T, class V>
class external_lin_writer : public external_writer_base<T, V> {
public:    
    external_lin_writer(const Rcpp::RObject& incoming) : external_writer_base<T, V>(incoming) {
        auto data_type=this->get_type();
        const char* type=data_type.c_str();
        auto classinfo=get_class_package(this->original);
        const char* pkg=classinfo.second.c_str();

        // Getting all required functions from the corresponding shared library.
        auto store_col2int_name=Set_external_name(classinfo, type, "output", "setCol", "integer");
        store_col_int=reinterpret_cast<void (*)(void *, size_t, RcppIntIt*, size_t, size_t)>(R_GetCCallable(pkg, store_col2int_name.c_str()));

        auto store_row2int_name=Set_external_name(classinfo, type, "output", "setRow", "integer");
        store_row_int=reinterpret_cast<void (*)(void *, size_t, RcppIntIt*, size_t, size_t)>(R_GetCCallable(pkg, store_row2int_name.c_str()));

        auto store_col2dbl_name=Set_external_name(classinfo, type, "output", "setCol", "numeric");
        store_col_dbl=reinterpret_cast<void (*)(void *, size_t, RcppNumIt*, size_t, size_t)>(R_GetCCallable(pkg, store_col2dbl_name.c_str()));

        auto store_row2dbl_name=Set_external_name(classinfo, type, "output", "setRow", "numeric");
        store_row_dbl=reinterpret_cast<void (*)(void *, size_t, RcppNumIt*, size_t, size_t)>(R_GetCCallable(pkg, store_row2dbl_name.c_str()));

        auto store_col2int_indexed_name=Set_external_name(classinfo, type, "output", "setColIndexed", "integer");
        store_col_indexed_int=reinterpret_cast<void (*)(void *, size_t, size_t, RcppIntIt*, RcppIntIt*)>(R_GetCCallable(pkg, store_col2int_name.c_str()));

        auto store_row2int_name=Set_external_name(classinfo, type, "output", "setRowIndexed", "integer");
        store_row_indexed_int=reinterpret_cast<void (*)(void *, size_t, size_t, RcppIntIt*, RcppIntIt*)>(R_GetCCallable(pkg, store_row2int_name.c_str()));

        auto store_col2dbl_name=Set_external_name(classinfo, type, "output", "setColIndexed", "numeric");
        store_col_indexed_dbl=reinterpret_cast<void (*)(void *, size_t, size_t, RcppIntIt*, RcppNumIt*)>(R_GetCCallable(pkg, store_col2dbl_name.c_str()));

        auto store_row2dbl_name=Set_external_name(classinfo, type, "output", "setRowIndexed", "numeric");
        store_row_indexed_dbl=reinterpret_cast<void (*)(void *, size_t, size_t, RcppIntIt*, RcppNumIt*)>(R_GetCCallable(pkg, store_row2dbl_name.c_str()));

        auto load_col2int_name=get_external_name(classinfo, type, "output", "getCol", "integer");
        load_col_int=reinterpret_cast<void (*)(void *, size_t, RcppIntIt*, size_t, size_t)>(R_GetCCallable(pkg, load_col2int_name.c_str()));

        auto load_row2int_name=get_external_name(classinfo, type, "output", "getRow", "integer");
        load_row_int=reinterpret_cast<void (*)(void *, size_t, RcppIntIt*, size_t, size_t)>(R_GetCCallable(pkg, load_row2int_name.c_str()));

        auto load_col2dbl_name=get_external_name(classinfo, type, "output", "getCol", "numeric");
        load_col_dbl=reinterpret_cast<void (*)(void *, size_t, RcppNumIt*, size_t, size_t)>(R_GetCCallable(pkg, load_col2dbl_name.c_str()));

        auto load_row2dbl_name=get_external_name(classinfo, type, "output", "getRow", "numeric");
        load_row_dbl=reinterpret_cast<void (*)(void *, size_t, RcppNumIt*, size_t, size_t)>(R_GetCCallable(pkg, load_row2dbl_name.c_str()));

        return;
    }

    ~external_lin_writer() = default;
    external_lin_writer(const external_lin_writer&) = default;
    external_lin_writer& operator=(const external_lin_writer&) = default;
    external_lin_writer(external_lin_writer&&) = default;
    external_lin_writer& operator=(external_lin_writer&&) = default;

    // Basic setters
    void set_row(size_t r, RcppIntIt out, size_t first, size_t last) {
        this->check_rowargs(r, first, last);
        store_row_int(this->ex.set(), r, &out, first, last);
        return;
    }

    void set_row(size_t r, RcppNumIt out, size_t first, size_t last) {
        this->check_rowargs(r, first, last);
        store_row_dbl(this->ex.set(), r, &out, first, last);
        return;
    }

    void set_col(size_t c, RcppIntIt out, size_t first, size_t last) {
        this->check_colargs(c, first, last);
        store_col_int(this->ex.set(), c, &out, first, last);
        return;
    }

    void set_col(size_t c, RcppNumIt out, size_t first, size_t last) {
        this->check_colargs(c, first, last);
        store_col_dbl(this->ex.set(), c, &out, first, last);
        return;
    }

    // Indexed setters
    void set_col_indexed(size_t c, size_t n, RcppIntIt idx, RcppIntIt in) {
        this->check_colargs(c);
        store_col_indexed_int(ex.get(), c, n, &idx, &in);
        return;
    }

    void set_col_indexed(size_t c, size_t n, RcppIntIt idx, RcppNumIt in) {
        this->check_colargs(c);
        store_col_indexed_dbl(ex.get(), c, n, &idx, &in);
        return;
    }

    void set_row_indexed(size_t c, size_t n, RcppIntIt idx, RcppIntIt in) {
        this->check_rowargs(c);
        store_row_indexed_int(ex.get(), c, n, &idx, &in);
        return;
    }
    
    void set_row_indexed(size_t c, size_t n, RcppIntIt idx, RcppNumIt in) {
        this->check_rowargs(c);
        store_row_indexed_dbl(ex.get(), c, n, &idx, &in);
        return;
    }

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


private:
    // Typedef'd for convenience.
    typedef Rcpp::IntegerVector::iterator RcppIntIt;
    typedef Rcpp::NumericVector::iterator RcppNumIt;

    void (*store_col_int) (void *, size_t, RcppIntIt*, size_t, size_t);
    void (*store_row_int) (void *, size_t, RcppIntIt*, size_t, size_t);
    void (*store_col_dbl) (void *, size_t, RcppNumIt*, size_t, size_t);
    void (*store_row_dbl) (void *, size_t, RcppNumIt*, size_t, size_t);

    void (*store_col_indexed_int) (void *, size_t, size_t, RcppIntIt*, RcppIntIt*);
    void (*store_row_indexed_int) (void *, size_t, size_t, RcppIntIt*, RcppIntIt*);
    void (*store_col_indexed_dbl) (void *, size_t, size_t, RcppIntIt*, RcppNumIt*);
    void (*store_row_indexed_dbl) (void *, size_t, size_t, RcppIntIt*, RcppNumIt*);

    void (*load_col_int) (void *, size_t, RcppIntIt*, size_t, size_t);
    void (*load_row_int) (void *, size_t, RcppIntIt*, size_t, size_t);
    void (*load_col_dbl) (void *, size_t, RcppNumIt*, size_t, size_t);
    void (*load_row_dbl) (void *, size_t, RcppNumIt*, size_t, size_t);
};

}

#endif
