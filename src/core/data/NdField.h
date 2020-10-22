//
// Created by rja on 21/10/2020.
//

#ifndef M7_NDFIELD_H
#define M7_NDFIELD_H

#include "Field.h"
#include "src/core/nd/NdFormat.h"

template<size_t nind>
struct NdFieldX : FieldX {
    NdFormat<nind> m_format;

    NdFieldX(TableX *table, std::array<size_t, nind> shape, size_t element_size, std::string description) :
            FieldX(table, NdFormat<nind>(shape).nelement(), element_size, description),
            m_format(shape) {}

    template<typename ...Args>
    char *raw_ptr(const size_t &irow, Args...inds) const {
        return FieldX::raw_ptr(irow, m_format.flat(inds...));
    }

    template<typename ...Args>
    std::pair<char *, size_t> raw_view(const size_t &irow, Args...inds) {
        return {raw_ptr(irow, inds...), m_element_size};
    }
};



#endif //M7_NDFIELD_H
