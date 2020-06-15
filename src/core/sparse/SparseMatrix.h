//
// Created by rja on 14/06/2020.
//

#ifndef M7_SPARSEMATRIX_H
#define M7_SPARSEMATRIX_H

#include <vector>

template <typename T>
struct SparseEntry {
    size_t icol;
    T element;
};

template<typename T>
class SparseMatrix {

    std::vector<std::vector<SparseEntry<T>>> m_data;

    SparseMatrix(){}

    void set(const size_t& irow, const size_t& icol, const T& element){
        if (irow>=m_data.size()) m_data.resize(irow+1);
        m_data[irow].push_back(SparseEntry<T>{icol, element});
    }

    //void apply(SparseVector)

};


#endif //M7_SPARSEMATRIX_H
