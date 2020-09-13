#ifndef BEACHMAT_CSPARSE_MATRIX_H
#define BEACHMAT_CSPARSE_MATRIX_H

#include "Rcpp.h"

#include "any_matrix.h"
#include "Csparse_core.h"
#include "utils.h"

#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

namespace beachmat {

template <typename T>
class Csparse_matrix : public any_matrix<T> {
public:
    Csparse_matrix () {}
    virtual ~Csparse_matrix() = default;
    Csparse_matrix(const Csparse_matrix&) = default;
    Csparse_matrix& operator=(const Csparse_matrix&) = default;
    Csparse_matrix(Csparse_matrix&&) = default;
    Csparse_matrix& operator=(Csparse_matrix&&) = default;

    virtual sparse_index<T, int> get_col(size_t, size_t, size_t);

    T* get_col(size_t r, T* work, size_t first, size_t last) {
        auto index = this->get_col(r, first, last);
        std::fill(work, work + last - first, static_cast<T>(0));
        for (size_t i=0; i < index.n; ++i) {
            work[index.i[i]] = index.x[i];
        }
        return work;
    }

    virtual sparse_index<T, int> get_row(size_t, T*, int*, size_t, size_t);
};

template <class V, typename T = typename V::stored_type>
class gCMatrix : public Csparse_matrix<T> {
public:
    virtual ~gCMatrix() = default;
    gCMatrix(const gCMatrix&) = default;
    gCMatrix& operator=(const gCMatrix&) = default;
    gCMatrix(gCMatrix&&) = default;
    gCMatrix& operator=(gCMatrix&&) = default;

    gCMatrix(Rcpp::RObject mat) { 
        this->fill_dims(get_safe_slot(mat, "Dim"));
        const size_t& NC=this->ncol;
        const size_t& NR=this->nrow;

        Rcpp::RObject temp_i=get_safe_slot(mat, "i");
        if (temp_i.sexp_type()!=INTSXP) { 
            auto ctype = get_class_name(mat);
            throw std::runtime_error(std::string("'i' slot in a ") + ctype + " object should be integer"); 
        }
        i=temp_i;

        Rcpp::RObject temp_p=get_safe_slot(mat, "p");
        if (temp_p.sexp_type()!=INTSXP) { 
            auto ctype = get_class_name(mat);
            throw std::runtime_error(std::string("'p' slot in a ") + ctype + " object should be integer");
        }
        p=temp_p;

        Rcpp::RObject temp_x=get_safe_slot(mat, "x");
        if (temp_x.sexp_type()!=x.sexp_type()) { 
            auto ctype = get_class_name(mat);
            throw std::runtime_error(std::string("'x' slot in a ") + ctype + " object should be " + translate_type(x.sexp_type()));
        }
        x=temp_x;

        if (x.size()!=i.size()) { 
            auto ctype = get_class_name(mat);
            throw std::runtime_error(std::string("'x' and 'i' slots in a ") + ctype + " object should have the same length"); 
        }
        if (NC+1!=static_cast<size_t>(p.size())) { 
            auto ctype = get_class_name(mat);
            throw std::runtime_error(std::string("length of 'p' slot in a ") + ctype + " object should be equal to 'ncol+1'"); 
        }
        if (p[0]!=0) { 
            auto ctype = get_class_name(mat);
            throw std::runtime_error(std::string("first element of 'p' in a ") + ctype + " object should be 0"); 
        }
        if (p[NC]!=x.size()) { 
            auto ctype = get_class_name(mat);
            throw std::runtime_error(std::string("last element of 'p' in a ") + ctype + " object should be 'length(x)'"); 
        }

        // Checking all the indices.
        auto pIt=p.begin();
        for (size_t px=0; px<NC; ++px) {
            const int& current=*pIt;
            if (current < 0) { 
                auto ctype = get_class_name(mat);
                throw std::runtime_error(std::string("'p' slot in a ") + ctype + " object should contain non-negative values"); 
            }
            if (current > *(++pIt)) { 
                auto ctype = get_class_name(mat);
                throw std::runtime_error(std::string("'p' slot in a ") + ctype + " object should be sorted"); 
            }
        }

        pIt=p.begin();
        for (size_t px=0; px<NC; ++px) {
            int left=*pIt; // Integers as that's R's storage type. 
            int right=*(++pIt)-1; // Not checking the last element, as this is the start of the next column.
            auto iIt=i.begin()+left;

            for (int ix=left; ix<right; ++ix) {
                const int& current=*iIt;
                if (current > *(++iIt)) {
                    auto ctype = get_class_name(mat);
                    throw std::runtime_error(std::string("'i' in each column of a ") + ctype + " object should be sorted");
                }
            }
        }

        for (auto iIt=i.begin(); iIt!=i.end(); ++iIt) {
            const int& curi=*iIt;
            if (curi<0 || static_cast<size_t>(curi)>=NR) {
                auto ctype = get_class_name(mat);
                throw std::runtime_error(std::string("'i' slot in a ") + ctype + " object should contain elements in [0, nrow)");
            }
        }

        core=Csparse_core<T, int, int>(i.size(), x.begin(), i.begin(), NR, NC, p.begin());
        return;                
    }

    sparse_index<T, int> get_col(size_t c, size_t first, size_t last) {
        return core.get_col(c, first, last);
    }

    sparse_index<T, int> get_row(size_t r, T* work_x, int* work_i, size_t first, size_t last) {
        return core.get_row(r, work_x, work_i, first, last);
    }

    T* get_row(size_t r, T* work, size_t first, size_t last) {
        core.get_row(r, work, first, last, 0);
        return work;
    }
private:
    Rcpp::IntegerVector i, p;
    V x;
    Csparse_core<T, int, int> core;
};

template <typename T>
struct sparse_triplet {
    sparse_triplet(int _i, int _j, T _x) : i(_i), j(_j), x(_x) {}
    int i;
    int j;
    T x;
};

template <class V, typename T = typename V::stored_type>
class sparse_seed : public Csparse_matrix<T> {
public:
    virtual ~sparse_seed() = default;
    sparse_seed(const sparse_seed&) = default;
    sparse_seed& operator=(const sparse_seed&) = default;
    sparse_seed(sparse_seed&&) = default;
    sparse_seed& operator=(sparse_seed&&) = default;

    sparse_seed(Rcpp::RObject seed) {
        this->fill_dims(get_safe_slot(seed, "dim"));
        const size_t& NC=this->ncol;
        const size_t& NR=this->nrow;
        p.resize(NC + 1);
    
        Rcpp::RObject temp_x=get_safe_slot(seed, "nzdata");
        if (temp_x.sexp_type()!=x.sexp_type()) { 
            auto ctype = get_class_name(seed);
            throw std::runtime_error(std::string("'nzdata' slot in a ") + ctype + " object should be " 
                + translate_type(x.sexp_type()));
        }
        x=temp_x;

        Rcpp::RObject temp_i=get_safe_slot(seed, "nzindex");
        if (temp_i.sexp_type() != INTSXP) { 
            auto ctype = get_class_name(seed);
            throw std::runtime_error(std::string("'nzindex' slot in a ") + ctype + " object should be integer"); 
        }
        Rcpp::IntegerMatrix temp_i2(temp_i);
        if (temp_i2.ncol() != 2) {
            auto ctype = get_class_name(seed);
            throw std::runtime_error(std::string("'nzindex' slot in a ") + ctype + " object should have two columns"); 
        }
        const size_t nnz = temp_i2.nrow();
        if (nnz != x.size()) {
            auto ctype = get_class_name(seed);
            throw std::runtime_error(std::string("incompatible 'nzindex' and 'nzdata' lengths in a ") + ctype + " object"); 
        }

        if (!nnz) {
            auto row_indices=temp_i2.column(0);
            auto col_indices=temp_i2.column(1);
            bool okay=true;
            {
                auto rowIt=row_indices.begin();
                auto colIt=col_indices.begin();

                for (size_t v = 0; v < nnz; ++v) {
                    auto lastR=*rowIt;
                    auto lastC=*colIt;
                    if (lastR <= 0 || lastR > NR || lastC <= 0 || lastC > NC) {
                        auto ctype = get_class_name(seed);
                        throw std::runtime_error(std::string("'nzindex' out of bounds in a ") + ctype + " object");
                    }

                    if (okay && v < nnz - 1) {
                        auto nextR=*(++rowIt);
                        auto nextC=*(++colIt);
                        if (lastC > nextC || (lastC == nextC && lastR > nextR)) {
                            okay = false;
                        }
                    }
                }
            }

            if (okay) {
                i = Rcpp::IntegerVector(row_indices.begin(), row_indices.end());
                for (auto& subi : i) { --subi; }

                auto colIt = col_indices.begin();
                auto colStart = colIt;
                for (int c = 0; c <= static_cast<int>(NC); ++c) {
                    while (colIt != col_indices.end() && *colIt <= c + 1) {
                        ++colIt;
                    }
                    p[c] = colIt - colStart;
                }
            } else {
                std::vector<sparse_triplet<T> > everything;
                everything.reserve(nnz);

                auto rowIt=row_indices.begin();
                auto colIt=col_indices.begin();
                auto xIt=x.begin();

                for (int v = 0; v < nnz; ++v, ++rowIt, ++colIt, ++xIt) {
                    everything.push_back(sparse_triplet<T>(*rowIt, *colIt, *xIt));
                }

                std::sort(everything.begin(), everything.end(), 
                    [] (const sparse_triplet<T>& left, const sparse_triplet<T>& right) -> bool {
                    if (left.j < right.j) {
                        return true;
                    } else if (left.j > right.j) {
                        return false;
                    } else if (left.i < right.i) {
                        return true;
                    } else if (left.i > right.i) {
                        return false;
                    } else if (left.x < right.x) {
                        return true;
                    } else {
                        return false;
                    }
                });

                i = Rcpp::IntegerVector(nnz);
                x = V(nnz);
                for (size_t v = 0; v < nnz; ++v) {
                    i[v] = everything[v].i - 1;
                    x[v] = everything[v].x;
                }

                auto eIt = everything.begin();
                auto eStart = eIt;
                for (int c = 0; c <= static_cast<int>(NC); ++c) {
                    while (eIt != everything.end() && eIt->j <= c + 1) {
                        ++eIt;
                    }
                    p[c] = eIt - eStart;
                }
            }
        }

        core=Csparse_core<T, int, size_t>(nnz, x.begin(), i.begin(), NR, NC, p.data());
        return;
    }

    sparse_index<T, int> get_col(size_t c, size_t first, size_t last) {
        return core.get_col(c, first, last);
    }

    sparse_index<T, int> get_row(size_t r, T* work_x, int* work_i, size_t first, size_t last) {
        return core.get_row(r, work_x, work_i, first, last);
    }

    T* get_row(size_t r, T* work, size_t first, size_t last) {
        core.get_row(r, work, first, last, sparse_seed<V, T>::get_empty());
        return work;
    }
private:
    Rcpp::IntegerVector i;
    V x;
    std::vector<size_t> p;
    Csparse_core<T, int, size_t> core;

    static const T get_empty();
};

template<>
const int sparse_seed<Rcpp::IntegerVector>::get_empty() { 
    return 0;
}

template<>
const int sparse_seed<Rcpp::LogicalVector>::get_empty() { 
    return 0;
}

template<>
const double sparse_seed<Rcpp::NumericVector>::get_empty() { 
    return 0;
}

template<>
const Rcpp::String sparse_seed<Rcpp::StringVector, Rcpp::String>::get_empty() {
    return Rcpp::String(0);
}

}

#endif
