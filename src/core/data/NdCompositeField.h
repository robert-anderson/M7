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

//    template<typename ...Args>
//    const char *raw_ptr(const size_t &icomponent, const size_t &irow, Args... inds) const {
//        return CompositeField::raw_ptr(irow, m_format.flatten(inds...));
//    }
//
//    template<typename ...Args>
//    std::pair<const char *, size_t> raw_view(const size_t &icomponent, const size_t &irow, Args... inds) const {
//        return CompositeField::raw_ptr(icomponent, irow, m_format.flatten(inds...));
//    }

//    template<typename ...Args>
//    defs::hash_t hash(const size_t &irow, Args... inds) const{
//        defs::hash_t res = 0;
//        for (size_t i=0ul; i<nfield(); ++i)
//            res^=m_fields[0]->hash(irow, m_format.flatten(inds...));
//        return res;
//    };

};

#endif //M7_NDCOMPOSITEFIELD_H
