#ifndef BEACHMAT_CSPARSE_MATRIX_H
#define BEACHMAT_CSPARSE_MATRIX_H

#include "Rcpp.h"

#include "any_matrix.h"
#include "utils.h"

#include <stdexcept>
#include <string>
#include <vector>

template <typename T>
class Csparse_matrix : any_matrix<T> {
public:
    virtual ~sparse_matrix() = default;
    sparse_matrix(const sparse_matrix&) = default;
    sparse_matrix& operator=(const sparse_matrix&) = default;
    sparse_matrix(sparse_matrix&&) = default;
    sparse_matrix& operator=(sparse_matrix&&) = default;

    T* get_col(size_t r, T* work, size_t first, size_t last) {
        auto index = get_col(r, first, last);
        std::fill(work, work + last - first, static_cast<T>(0));
        for (size_t i=0; i < index.n; ++i) {
            work[index.i[i]] = index.x[i];
        }
        return work;
    }

    virtual sparse_index<T, int> get_col(size_t, size_t, size_t);
    virtual sparse_index<T, int> get_col(size_t, sparse_index<T, int>&, size_t, size_t);
};

template <typename T, class V>
class gCMatrix : Csparse_matrix<T> {
public:
    virtual ~gcMatrix() = default;
    gcMatrix(const gcMatrix&) = default;
    gcMatrix& operator=(const gcMatrix&) = default;
    gcMatrix(gcMatrix&&) = default;
    gcMatrix& operator=(gcMatrix&&) = default;

    gcMatrix(Rcpp::RObject mat) { 
        this->fill_dims(get_safe_slot(incoming, "Dim"));
        const size_t& NC=this->ncol;
        const size_t& NR=this->nrow;

        Rcpp::RObject temp_i=get_safe_slot(incoming, "i");
        if (temp_i.sexp_type()!=INTSXP) { 
            throw std::runtime_error(std::string("'i' slot in a ") + ctype + " object should be integer"); 
        }
        i=temp_i;

        Rcpp::RObject temp_p=get_safe_slot(incoming, "p");
        if (temp_p.sexp_type()!=INTSXP) { 
            throw std::runtime_error(std::string("'p' slot in a ") + ctype + " object should be integer");
        }
        p=temp_p;

        Rcpp::RObject temp_x=get_safe_slot(incoming, "x");
        if (temp_x.sexp_type()!=x.sexp_type()) { 
            throw std::runtime_error(std::string("'x' slot in a ") + ctype + " object should be " + translate_type(x.sexp_type()));
        }
        x=temp_x;

        if (x.size()!=i.size()) { 
            throw std::runtime_error(std::string("'x' and 'i' slots in a ") + ctype + " object should have the same length"); 
        }
        if (NC+1!=static_cast<size_t>(p.size())) { 
            throw std::runtime_error(std::string("length of 'p' slot in a ") + ctype + " object should be equal to 'ncol+1'"); 
        }
        if (p[0]!=0) { 
            throw std::runtime_error(std::string("first element of 'p' in a ") + ctype + " object should be 0"); 
        }
        if (p[NC]!=x.size()) { 
            throw std::runtime_error(std::string("last element of 'p' in a ") + ctype + " object should be 'length(x)'"); 
        }

        // Checking all the indices.
        auto pIt=p.begin();
        for (size_t px=0; px<NC; ++px) {
            const int& current=*pIt;
            if (current < 0) { 
                throw std::runtime_error(std::string("'p' slot in a ") + ctype + " object should contain non-negative values"); 
            }
            if (current > *(++pIt)) { 
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
                    throw std::runtime_error(std::string("'i' in each column of a ") + ctype + " object should be sorted");
                }
            }
        }

        for (auto iIt=i.begin(); iIt!=i.end(); ++iIt) {
            const int& curi=*iIt;
            if (curi<0 || static_cast<size_t>(curi)>=NR) {
                throw std::runtime_error(std::string("'i' slot in a ") + ctype + " object should contain elements in [0, nrow)");
            }
        }

        core=Sparse_core<T, int, int>(i.size(), x.begin(), i.begin(), NR, NC, p.begin());
        return;                
    }

    sparse_index<T, int> get_col(size_t c, size_t first, size_t last) {
        return core.get_col(c, first, last);
    }

    sparse_index<T, int> get_row(size_t r, sparse_index<T, int>& work, size_t first, size_t last) {
        return core.get_row(r, work, first, last);
    }

    void get_row(size_t r, T* work, size_t first, size_t last) {
        return core.get_row(r, work, first, last);
    }
private:
    Rcpp::IntegerVector i, p;
    V x;
    Csparse_core<T, int, int> core;
};

template <typename T, class V>
class sparse_seed : Csparse_matrix<T> {
private:
    struct sparse_triplet {
        sparse_triplet(int _i, int _j, T _x) : i(_i), j(_j), x(_x) {}
        int i;
        int j;
        T x;
    };
public:
    virtual ~sparse_seed() = default;
    sparse_seed(const sparse_seed&) = default;
    sparse_seed& operator=(const sparse_seed&) = default;
    sparse_seed(sparse_seed&&) = default;
    sparse_seed& operator=(sparse_seed&&) = default;

    sparse_seed(Rcpp::RObject seed) {
        this->fill_dims(get_safe_slot(incoming, "dim"));
        const size_t& NC=this->ncol;
        const size_t& NR=this->nrow;
        p.resize(NC + 1);
    
        Rcpp::RObject temp_x=get_safe_slot(incoming, "nzdata");
        if (temp_x.sexp_type()!=x.sexp_type()) { 
            throw std::runtime_error(std::string("'nzdata' slot in a ") + ctype + " object should be " 
                + translate_type(x.sexp_type()));
        }
        x=temp_x;

        Rcpp::RObject temp_i=get_safe_slot(incoming, "nzindex");
        if (temp_i.sexp_type() != INTSXP) { 
            throw std::runtime_error(std::string("'nzindex' slot in a ") + ctype + " object should be integer"); 
        }
        Rcpp::IntegerMatrix temp_i2 = temp_i;
        if (temp_i2.ncol() != 2) {
            throw std::runtime_error(std::string("'nzindex' slot in a ") + ctype + " object should have two columns"); 
        }
        const size_t nnz = temp_i2.nrow();
        if (nnz != x.size()) {
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
                        throw std::runtime_error(std::string("'nzindex' out of bounds in a ") + ctype + " object");
                    }

                    if (okay && r < nnz - 1) {
                        auto nextR=*(++rowIt);
                        auto nextR=*(++colIt);
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
                std::vector<sparse_triplet> everything;
                everything.reserve(nnz);

                auto rowIt=row_indices.begin();
                auto colIt=col_indices.begin();
                auto xIt=x.begin();

                for (int v = 0; v < nnz; ++v, ++rowIt, ++colIt, ++xIt) {
                    everything.push_back(sparse_triplet(*rowIt, *colIt, *xIt));
                }

                std::sort(everything.begin(), everything.end(), 
                    [] (const sparse_triplet& left, const sparse_triplet& right) -> bool {
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
                })

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

        core=Csparse_core<T, int, size_t>(nnz, i.begin(), x.begin(), NR, NC, p.data());
    }
private:
    Rcpp::IntegerVector i;
    V x;
    std::vector<size_t> p;
    Csparse_core<T, int, size_t> core;
}

#endif
