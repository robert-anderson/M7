//
// Created by RJA on 18/09/2020.
//

#ifndef M7_ARRAYFORMAT_H
#define M7_ARRAYFORMAT_H

#include "Indexer_.h"

template<size_t nind>
struct ArrayFormat : ArrayFormatBase<nind> {
    template<typename ...Args>
    ArrayFormat(const size_t &first, Args ...shape) : ArrayFormatBase<nind>(first, shape...) {}

    Indexer_<nind, 1> operator[](size_t i) {
        return Indexer_<nind, 1>(*this, i);
    }
};



#endif //M7_ARRAYFORMAT_H
