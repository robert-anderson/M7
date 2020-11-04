//
// Created by RJA on 28/10/2020.
//

#ifndef M7_FIELDS_H
#define M7_FIELDS_H

#include "NumericSpecifier.h"
#include "NumericArraySpecifier.h"
#include "BitsetSpecifier.h"
#include "DeterminantSpecifier.h"
#include "BosonOnvSpecifier.h"
#include "FermionBosonState.h"
#include "Flag.h"

/*
 * some convenient wrappers for the field types
 */

namespace fields {

    template<size_t nind>
    using Flags = NdFlag<nind>;

    using Flag = Flags<0ul>;

    template<typename T, size_t nind>
    struct Numbers: NdField<NumericSpecifier<T>, nind>{
        template<typename ...Args>
        Numbers(TableX *table, std::string description, Args... shape) :
            NdField<NumericSpecifier<T>, nind>(table, {}, description, shape...) {}
    };

    template<typename T>
    using Number = Numbers<T, 0ul>;

    template<typename T, size_t nind, size_t nind_element>
    using NumberArrays = NdField<NumericArraySpecifier<T, nind_element>, nind>;

    template<size_t nind>
    using Bitsets = NdField<BitsetSpecifier, nind>;

    using Bitset = Bitsets<0ul>;

    template<size_t nind>
    using Determinants = NdField<DeterminantSpecifier, nind>;

    using Determinant = Determinants<0ul>;

    template<size_t nind>
    using BosonOnvs = NdField<BosonOnvSpecifier, nind>;

    using BosonOnv = BosonOnvs<0ul>;

    template<size_t nind>
    using FermionBosonConfigurations = fb_state::Field<nind>;

    using FermionBosonConfiguration = FermionBosonConfigurations<0ul>;

    template<size_t nind, bool bosons>
    struct ConfigurationSelector {};

    template<size_t nind>
    struct ConfigurationSelector<nind, false> {
        typedef Determinants<nind> type;
    };

    template<size_t nind>
    struct ConfigurationSelector<nind, true> {
        typedef FermionBosonConfigurations<nind> type;
    };

    template<size_t nind>
    using Configurations = typename ConfigurationSelector<nind, defs::bosons>::type;

    using Configuration = Configurations<0ul>;
}

#endif //M7_FIELDS_H
