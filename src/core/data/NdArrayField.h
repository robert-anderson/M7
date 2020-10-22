//
// Created by rja on 21/10/2020.
//

#ifndef M7_NDARRAYFIELD_H
#define M7_NDARRAYFIELD_H

#include "NdField.h"

template<size_t nind, size_t nind_view>
struct NdArrayField : NdFieldX<nind> {
    static_assert(nind_view, "If the view is scalar, use NdField instead");
    NdFormat<nind_view> m_view_format;

    NdArrayField(TableX *table, std::array<size_t, nind> shape,
                 std::array<size_t, nind_view> view_shape, size_t element_size, std::string description) :
            NdFieldX<nind>(table, shape, NdFormat<nind_view>(view_shape).nelement() * element_size, description),
            m_view_format(view_shape) {}

};



#endif //M7_NDARRAYFIELD_H
