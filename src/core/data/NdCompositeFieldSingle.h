//
// Created by rja on 01/11/2020.
//

#ifndef M7_NDCOMPOSITEFIELDSINGLE_H
#define M7_NDCOMPOSITEFIELDSINGLE_H

#include "NdCompositeField.h"
#include "NdField.h"
#include "Table.h"

template<typename field_t, size_t nind>
struct NdCompositeFieldSingle : NdCompositeField<nind> {
    NdFieldX<field_t, nind> m_field;

    NdCompositeFieldSingle(TableX *table, field_t&& field, std::string description, NdFormat<nind> format) :
            NdCompositeField<nind>(table, format),
            m_field(table, std::move(field), description, format) {}

    typedef typename field_t::view_t view_t;
    typedef typename field_t::const_view_t const_view_t;

    template<typename ...Args>
    view_t operator()(const size_t& irow, Args... inds){
        return m_field(irow, inds...);
    }

    template<typename ...Args>
    const_view_t operator()(const size_t& irow, Args... inds) const{
        return m_field(irow, inds...);
    }
};



#endif //M7_NDCOMPOSITEFIELDSINGLE_H
