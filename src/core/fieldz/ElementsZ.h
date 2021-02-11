//
// Created by rja on 10/02/2021.
//

#ifndef M7_ELEMENTSZ_H
#define M7_ELEMENTSZ_H

#include "NdMultiFieldZ.h"
#include "BufferedTableZ.h"

/*
template<size_t nind, typename ...Args>
struct SingleFieldRowZ : RowZ {
    NdMultiFieldZ<nind, ...Args> m_field;
    SingleFieldRowZ(std::array<size_t, nind> shape, Args&&... subfields) :
    RowZ(), m_field(this, shape, subfields...) {}
};

template<typename nd_field_t>
struct BufferedSingleFieldTableZ : BufferedTableZ<SingleFieldRowZ<nd_field_t>> {
    template<typename ...Args>
    SingleFieldTableZ(Args... args):
            TableZ<SingleFieldRowZ<nd_field_t>>(SingleFieldRowZ<nd_field_t>(args...)){}
};

template<typename nd_field_t>
struct ElementZ : nd_field_t, BufferedTableZ<row_t> {
    ElementZ(Args... args): NdFieldZ<0ul, row_t>()
};
*/

#endif //M7_ELEMENTSZ_H
