#ifndef BEACHMAT_CSPARSE_CORE_H
#define BEACHMAT_CSPARSE_CORE_H

#include "Rcpp.h"
#include <algorithm>
#include <vector>

namespace beachmat {

template <typename T, typename I>
struct sparse_index {
    sparse_index(size_t _n, const T* _x, const I* _i) : n(_n), x(_x), i(_i) {}
    size_t n;
    const T* x;
    const I* i;
};

template <typename T, typename I, typename P>
class Csparse_core {
public:
    Csparse_core();
    Csparse_core(const size_t _n, const T* _x, const I* _i, const size_t _nr, const size_t _nc, const P* _p) : 
        n(_n), nr(_nr), nc(_nc), x(_x), i(_i), p(_p), currow(0), curstart(0), curend(nc) {}
    
    sparse_index get_col(size_t c, size_t first, size_t last) {
        check_colargs(c, first, last);
        const auto pstart=p[c]; 
        auto iIt=i.begin()+pstart, 
             eIt=i.begin()+p[c+1]; 
        auto xIt=x.begin()+pstart;

        if (first) { // Jumping ahead if non-zero.
            auto new_iIt=std::lower_bound(iIt, eIt, first);
            xIt+=(new_iIt-iIt);
            iIt=new_iIt;
        } 

        if (last!=nr) { // Jumping to last element.
            eIt=std::lower_bound(iIt, eIt, last);
        }

        return sparse_index<T, I>(eIt - iIt, xIt, iIt); 
    }

    void get_row(size_t r, T* work, size_t first, size_t last) {
        check_rowargs(r, first, last);
        update_indices(r, first, last);
        std::fill(work, work + last - first, static_cast<T>(0));

        auto pIt = p + first + 1; // Points to first-past-the-end for each 'c'.
        for (size_t c = first; c < last; ++c, ++pIt, ++work) { 
            const int& idex = indices[c];
            if (idex != *pIt && static_cast<size_t>(i[idex]) == r) { 
                (*work)=x[idex]; 
            }
        } 
        return;  
    }

    void get_row(size_t r, sparse_index& work, size_t first, size_t last) {
        check_rowargs(r, first, last);
        update_indices(r, first, last);

        auto pIt = p + first + 1; // Points to first-past-the-end for each 'c'.
        auto& counter = (work.n = 0);

        for (size_t c = first; c < last; ++c, ++pIt, ++work) { 
            const int& idex = indices[c];
            if (idex != *pIt && static_cast<size_t>(i[idex]) == r) { 
                work.i[counter] = i[idex];
                work.x[counter] = x[idex];
                ++counter;
            }
        }
        return;  
    }
private:
    const size_t n, nr, nc;
    typename const T* x;
    typename const I* i;
    typename const P* p;

    size_t currow, curstart, curend;
    std::vector<P> indices; 

    void update_indices(size_t r, size_t first, size_t last) {
        /* Initializing the indices upon the first request, assuming currow=0 based on initialization above.
         * This avoids using up space for the indices if we never do row access.
         */
        if (indices.size() != nr) {
            indices = std::vector<int>(p, p + nc);
        }

        /* If left/right slice are not equal to what is stored, we reset the indices,
         * so that the code below will know to recompute them. It's too much effort
         * to try to figure out exactly which columns need recomputing; just do them all.
         */
        if (first != curstart || last != curend) {
            std::copy(p, p + nc, indices.begin());
            currow=0;
        }

        /* entry of 'indices' for each column should contain the index of the first
         * element with row number not less than 'r'. If no such element exists, it
         * will contain the index of the first element of the next column.
         */
        if (r == currow) { 
            return; 
        } 

        const P* pIt = p + first;
        if (r == currow+1) {
            ++pIt; // points to the first-past-the-end element, at any given 'c'.
            for (size_t c=first; c<last; ++c, ++pIt) {
                P& curdex = indices[c];
                if (curdex != *pIt && static_cast<size_t>(i[curdex]) < r) { 
                    ++curdex;
                }
            }
        } else if (r+1 == currow) {
            for (size_t c=first; c<last; ++c, ++pIt) {
                P& curdex = indices[c];
                if (curdex != *pIt && static_cast<size_t>(i[curdex-1]) >= r) { 
                    --curdex;
                }
            }

        } else { 
            if (r > currow) {
                ++pIt; // points to the first-past-the-end element, at any given 'c'.
                for (size_t c = first; c < last; ++c, ++pIt) { 
                    indices[c] = std::lower_bound(i + curdex, i + *pIt, r) - i;
                }
            } else { 
                for (size_t c = first; c < last; ++c, ++pIt) {
                    indices[c] = std::lower_bound(i + *pIt, i + indices[c], r) - i;
                }
            }
        }

        currow=r;
        return;
    }
};

}

#endif
