//
// Created by rja on 14/06/2020.
//

#ifndef M7_SPARSEMATRIX_H
#define M7_SPARSEMATRIX_H

#include <vector>
#include <src/core/thread/Atomic.h>
#include "src/core/linalg/Matrix.h"
#include "omp.h"

template <typename T>
struct SparseEntry {
    size_t icol;
    T element;
};

template<typename T>
class SparseMatrix {
    std::vector<std::vector<SparseEntry<T>>> m_data;

public:
    SparseMatrix(){}

    size_t nrow(){return m_data.size();}

    void resize(const size_t nrow){
        ASSERT(nrow>m_data.size());
        m_data.resize(nrow);
    }

    void expand(const size_t delta_nrow){
        resize(m_data.size()+delta_nrow);
    }

    T& operator()(const size_t& irow, const size_t& icol){
        if (irow>=m_data.size()) {
            ASSERT(!omp_get_level());
            resize(irow+1);
        }
        m_data[irow].push_back(SparseEntry<T>{icol, 0});
        return m_data[irow].back().element;
    }

    bool empty(){return m_data.empty();}

    void multiply(const std::vector<T> &in, std::vector<T> &out){
#pragma omp parallel for default(none) shared(in, out)
        for (size_t irow=0; irow<m_data.size(); ++irow){
            for (auto entry=m_data[irow].begin(); entry!=m_data[irow].end(); entry++){
                as_atomic(out[entry->icol])+=entry->element*in[irow];
            }
        }
    }

};


#endif //M7_SPARSEMATRIX_H
