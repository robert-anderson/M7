//
// Created by rja on 21/10/2020.
//

#ifndef M7_NDARRAYFIELD_H
#define M7_NDARRAYFIELD_H

#include "Field.h"

#if 0
template<size_t nind, size_t nind_view>
struct NdArrayField : NdFieldX<nind> {
    static_assert(nind_view, "If the view is scalar, use NdField instead");
    NdFormat<nind_view> m_view_format;

    NdArrayField(TableX *table, std::array<size_t, nind> shape,
                 std::array<size_t, nind_view> view_shape, size_t element_size,
                 const std::type_info& type_info, std::string description) :
            NdFieldX<nind>(table, shape, NdFormat<nind_view>(view_shape).nelement() * element_size,
                           type_info, description),
            m_view_format(view_shape) {}

    struct View : NdFieldX<nind>::View {
        template<typename ...Args>
        View(const NdFieldX<nind>& field, const size_t &irow, const size_t& iflat):
                NdFieldX<nind>::View(field, irow, iflat) {}
    };

    std::map<std::string, std::string> details() const override {
        auto map = NdFieldX<nind>::details();
        map["element rank"] = std::to_string(nind_view);
        if (nind_view) map["element shape"] = utils::to_string(m_view_format.shape());
        return map;
    }
};


#endif //M7_NDARRAYFIELD_H
#endif //M7_NDARRAYFIELD_H
