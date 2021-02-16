//
// Created by rja on 16/02/2021.
//

#ifndef M7_NDFIELDBASEZ_H
#define M7_NDFIELDBASEZ_H

#include "FieldBaseZ.h"
#include "src/core/nd/NdFormat.h"



template<size_t nind_item>
struct ItemFormattedFieldBaseZ : FieldBaseZ {
    NdFormat<nind_item> *m_item_format = nullptr;
    ItemFormattedFieldBaseZ(size_t item_size, const std::type_info &type_info):
    FieldBaseZ(item_size, type_info){}
};

template<typename T, size_t nind_item, size_t nind_element>
struct FullyFormattedFieldBaseZ : ItemFormattedFieldBaseZ<nind_item> {
    using ItemFormattedFieldBaseZ<nind_item>::m_item_format;
    NdFormat<nind_element> m_element_format;
    FullyFormattedFieldBaseZ(size_t item_size, std::array<size_t, nind_element> element_shape):
            ItemFormattedFieldBaseZ<nind_item>(item_size, typeid(T)), m_element_format(element_shape){}

    size_t nelement_all() const {
        return m_item_format->nelement()*m_element_format.nelement();
    }
};

#endif //M7_NDFIELDBASEZ_H
