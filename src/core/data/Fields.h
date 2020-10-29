//
// Created by RJA on 28/10/2020.
//

#ifndef M7_FIELDS_H
#define M7_FIELDS_H

#include "NumericField.h"
#include "NumericArrayField.h"
#include "BitsetField.h"
#include "DeterminantField.h"

/*
 * some convenient wrappers for the field types
 */

namespace fields {
    template<typename T>
    struct Number : NdFieldX<NumericFieldX<T>, 0> {
        Number(TableX* table, std::string description):
        NdFieldX<NumericFieldX<T>, 0>(table, {}, description){}
    };

    template<typename T, size_t nind>
    struct Numbers : NdFieldX<NumericFieldX<T>, nind> {
        template <typename...Args>
        Numbers(TableX* table, std::string description, Args...shape):
                NdFieldX<NumericFieldX<T>, nind>(table, {}, description, shape...){}
    };

    template<typename T, size_t nind_view>
    using NumberArray = NdFieldX<NumericArrayField<T, nind_view>, 0>;

    template<typename T, size_t nind_view, size_t nind_field>
    using NumberArrays = NdFieldX<NumericArrayField<T, nind_view>, nind_field>;

    using Bitset = NdFieldX<BitsetFieldX, 0>;

    template<size_t nind>
    using Bitsets = NdFieldX<BitsetFieldX, nind>;

    using Determinant = NdFieldX<DeterminantFieldX, 0>;

    template<size_t nind>
    using Determinants = NdFieldX<DeterminantFieldX, nind>;

//    using BosonOnv = NdFieldX<BosonOnvX, 0>;

//    template<size_t nind>
//    using BosonOnvs = NdFieldX<BosonOnvX, nind>;

}


#endif //M7_FIELDS_H
