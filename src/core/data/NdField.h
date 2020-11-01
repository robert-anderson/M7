//
// Created by rja on 21/10/2020.
//

#ifndef M7_NDFIELD_H
#define M7_NDFIELD_H

#include <src/core/util/utils.h>
#include "NdFieldBase.h"
#include "src/core/nd/NdFormat.h"
#include "NdCompositeField.h"

template<typename field_t, size_t nind>
struct NdFieldX : NdFieldBaseX {
    static_assert(std::is_base_of<FieldBaseX, field_t>::value,
                  "Field parameter arg must be derived from FieldBase");
    field_t m_field;
    NdFormat<nind> m_format;

    std::string element_string(char *ptr) const override {
        return static_cast<const FieldBaseX &>(m_field).element_string(ptr);
    }

    NdFieldX(TableX *table, field_t &&field, std::string description, NdFormat<nind> format):
            NdFieldBaseX(table, std::move(field), format.nelement(), description),
            m_field(std::move(field)), m_format(format) {
        m_details["field dimensionality"] = std::to_string(nind);
        if (nind) m_details["field shape"] = utils::to_string(m_format.shape());
    }

    NdFieldX(NdCompositeField<nind> *composite, field_t &&field, std::string description):
            NdFieldX(composite->m_table, std::move(field), description, composite->m_format){}

    template<typename ...Args>
    typename field_t::view_t operator()(const size_t &irow, Args... inds) {
        return m_field(raw_ptr(irow, m_format.flatten(inds...)));
    }

    template<typename ...Args>
    typename field_t::const_view_t operator()(const size_t &irow, Args... inds) const {
        return m_field(raw_ptr(irow, m_format.flatten(inds...)));
    }
};

#endif //M7_NDFIELD_H
