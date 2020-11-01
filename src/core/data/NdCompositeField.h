//
// Created by RJA on 01/11/2020.
//

#ifndef M7_NDCOMPOSITEFIELD_H
#define M7_NDCOMPOSITEFIELD_H

#include "CompositeField.h"
#include "src/core/nd/NdFormat.h"


template<size_t nind>
struct NdCompositeField : CompositeField {

    NdFormat<nind> m_format;

    NdCompositeField(TableX *tableX, NdFormat<nind> format) :
            CompositeField(tableX), m_format(format) {}

};

#endif //M7_NDCOMPOSITEFIELD_H
