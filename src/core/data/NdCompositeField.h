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

    template<typename ...Args>
    const char *raw_ptr(const size_t &icomponent, const size_t &irow, Args... inds) const {
        return CompositeField::raw_ptr(irow, m_format.flatten(inds...));
    }

    template<typename ...Args>
    std::pair<const char *, size_t> raw_view(const size_t &icomponent, const size_t &irow, Args... inds) const {
        return CompositeField::raw_ptr(icomponent, irow, m_format.flatten(inds...));
    }

};

#endif //M7_NDCOMPOSITEFIELD_H
