//
// Created by rja on 21/10/2020.
//

#ifndef M7_NUMERICFIELD_H
#define M7_NUMERICFIELD_H

#include "NdField.h"

template<typename T, size_t nind>
struct NumericField : NdField<nind> {
    NumericField(Table *table, std::array<size_t, nind> shape, std::string description) :
            NdField<nind>(table, shape, sizeof(T), description) {}

    template<typename ...Args>
    T &operator()(const size_t &irow, Args... inds) {
        return (T *) raw_ptr(irow, inds...);
    }
};


#endif //M7_NUMERICFIELD_H
