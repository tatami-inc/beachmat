#ifndef AARON_MATRIX_H
#define AARON_MATRIX_H

#include "Rcpp.h"

/* NOTE: I have not checked for valid inputs here to keep this demonstration code simple.
 * Real-life applications should add bound checks, otherwise segmentation faults may occur.
 * The recommended option is to derive from beachmat::dim_checker and add check_args() calls to each get_* method. 
 * See beachmat/simple_reader.h for an example.
 */

template<typename T, class V, class M>
class AaronMatrix {
public:
    AaronMatrix(const Rcpp::RObject& incoming) : mat(Rcpp::RObject(Rcpp::S4(incoming).slot("data"))) {}

    int get_nrow() const { return mat.nrow(); }
    int get_ncol() const { return mat.ncol(); }

    // Basic getters.

    T get(int r, int c) { return mat(r, c); }

    template<class Iter>
    void get_row(int r, Iter out, int first, int last) {
        auto currow=mat.row(r);
        std::copy(currow.begin()+first, currow.begin()+last, out);
        return;
    }

    template<class Iter>
    void get_col(int c, Iter out, int first, int last) {
        auto curcol=mat.column(c);
        std::copy(curcol.begin()+first, curcol.begin()+last, out);
        return;
    }

    // Multi getters.

    template<class Iter>
    void get_rows(int * r, int n, Iter out, int first, int last) {
        for (size_t c=first; c<last; ++c) {
            auto curcol=mat.column(c);
            auto it=curcol.begin();
            auto r_copy=r;
            for (size_t i=0; i<n; ++i, ++out, ++r_copy) {  
                (*out)=*(it + *r_copy);
            }
        }
    }

    template<class Iter>
    void get_cols(int* c, int n, Iter out, int first, int last) {
        for (int i=0; i<n; ++i) {
            get_col(*c, out, first, last);
            out+=last - first;
            ++c;
        }
        return;
    }

private:
    M mat;
};

#endif
