//
// Created by RJA on 28/10/2020.
//

#ifndef M7_FIELDS_H
#define M7_FIELDS_H

#include "NdCompositeFieldSingle.h"
#include "NumericField.h"
#include "NumericArrayField.h"
#include "BosonOnvField.h"
#include "BitsetField.h"
#include "DeterminantField.h"
#include "ConfigurationField.h"

/*
 * some convenient wrappers for the field types
 */

namespace fields {

    /*
     * single composite fields...
     */

    template<typename T>
    struct Number : NdCompositeFieldSingle<NumericFieldX<T>, 0ul> {
        Number(TableX *table, std::string description) :
                NdCompositeFieldSingle<NumericFieldX<T>, 0ul>(table, {}, description, {}) {}
    };

    template<typename T, size_t nind>
    struct Numbers : NdCompositeFieldSingle<NumericFieldX<T>, nind> {
        template<typename...Args>
        Numbers(TableX *table, std::string description, Args...shape):
                NdCompositeFieldSingle<NumericFieldX<T>, nind>(table, {}, description, {shape...}) {}
    };

    template<typename T, size_t nind_view>
    using NumberArray = NdCompositeFieldSingle<NumericArrayField<T, nind_view>, 0ul>;

    template<typename T, size_t nind_view, size_t nind_field>
    using NumberArrays = NdCompositeFieldSingle<NumericArrayField<T, nind_view>, nind_field>;

    using Bitset = NdCompositeFieldSingle<BitsetFieldX, 0ul>;

    template<size_t nind>
    using Bitsets = NdCompositeFieldSingle<BitsetFieldX, nind>;

    using Determinant = NdCompositeFieldSingle<DeterminantFieldX, 0ul>;

    template<size_t nind>
    using Determinants = NdCompositeFieldSingle<DeterminantFieldX, nind>;

    using BosonOnv = NdCompositeFieldSingle<BosonOnvField, 0ul>;

    template<size_t nind>
    using BosonOnvs = NdCompositeFieldSingle<BosonOnvField, nind>;


    /*
     * true composite fields...
     */

    template<size_t nind, bool bosons>
    struct ConfigurationSelector {
    };

    template<size_t nind>
    struct ConfigurationSelector<nind, false> {
        typedef FermionField <nind> type;
    };

    template<size_t nind>
    struct ConfigurationSelector<nind, true> {
        typedef FermionBosonField <nind> type;
    };

    using Configuration = typename ConfigurationSelector<0ul, defs::bosons>::type;

    template<size_t nind>
    using Configurations = typename ConfigurationSelector<nind, defs::bosons>::type;

}


#endif //M7_FIELDS_H
