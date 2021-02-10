//
// Created by rja on 09/02/2021.
//

#ifndef M7_FIELDSZ_H
#define M7_FIELDSZ_H

#include "NumberFieldZ.h"
#include "BitsetFieldZ.h"
#include "NdMultiFieldZ.h"

namespace fieldsz {

    template<typename T, size_t nind, size_t nind_element>
    struct NumberArrays : NdMultiFieldZ<nind, NumberFieldZ<T, nind_element>> {
        NumberArrays(RowZ* row, std::array<size_t, nind> shape, std::array<size_t, nind_element> element_shape):
                NdMultiFieldZ<nind, NumberFieldZ<T, nind_element>>(row, shape, element_shape) {}
    };


}


#endif //M7_FIELDSZ_H
