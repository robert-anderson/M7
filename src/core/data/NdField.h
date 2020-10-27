//
// Created by rja on 21/10/2020.
//

#ifndef M7_NDFIELD_H
#define M7_NDFIELD_H

#include <src/core/util/utils.h>
#include "NdFieldBase.h"
#include "src/core/nd/NdFormat.h"

template<typename field_t, size_t nind>
struct NdFieldX : NdFieldBaseX {
    static_assert(std::is_base_of<FieldBaseX, field_t>::value,
                  "Field parameter arg must be derived from FieldBase");
    field_t m_field;
    NdFormat<nind> m_format;

    std::string element_string(char *ptr) const override {
        return static_cast<const FieldBaseX &>(m_field).element_string(ptr);
    }

    template<typename ...Args>
    NdFieldX(TableX *table, field_t &&field, std::string description, Args... shape):
            NdFieldBaseX(table, NdFormat<nind>(shape...).nelement(),
                         static_cast<FieldBaseX &&>(field).m_element_size, typeid(field_t), description),
            m_field(std::move(field)), m_format(shape...) {}

    template<typename ...Args>
    typename field_t::view_t operator()(const size_t &irow, Args... inds) {
        return m_field(raw_ptr(irow, m_format.flatten(inds...)));
    }

    template<typename ...Args>
    typename field_t::const_view_t operator()(const size_t &irow, Args... inds) const {
        return m_field(raw_ptr(irow, m_format.flatten(inds...)));
    }
};

//template<size_t nind>
//struct NdFieldX : FieldX {
//    NdFormat<nind> m_format;
//
//    NdFieldX(TableX *table, std::array<size_t, nind> shape, size_t element_size,
//          const std::type_info &type_info, std::string description) :
//            FieldX(table, NdFormat<nind>(shape).nelement(), element_size, type_info, description),
//            m_format(shape) {}
//
//    struct View : FieldX::View {
//        View(const NdFieldX &field, const size_t &irow, const size_t &iflat) :
//                FieldX::View(field, irow, iflat) {}
//    };
//
//    template<typename ...Args>
//    char *raw_ptr(const size_t &irow, Args...inds) const {
//        return FieldX::raw_ptr(irow, m_format.flatten(inds...));
//    }
//
//    template<typename ...Args>
//    std::pair<char *, size_t> raw_view(const size_t &irow, Args...inds) {
//        return {raw_ptr(irow, inds...), m_element_size};
//    }
//
//    std::map<std::string, std::string> details() const override {
//        auto map = FieldX::details();
//        map["field rank"] = std::to_string(nind);
//        if (nind) map["field shape"] = utils::to_string(m_format.shape());
//        return map;
//    }
//};


#endif //M7_NDFIELD_H
