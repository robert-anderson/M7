//
// Created by RJA on 28/10/2020.
//

#ifndef M7_FIELDS_H
#define M7_FIELDS_H

#include "Specifiers.h"
#include "FermionBosonOnv.h"
#include "FlagField.h"

/*
 * some convenient wrappers for the field types
 */

namespace fields {

    template <typename spec_t>
    using Flags = FlagField<spec_t>;

    template<typename T, size_t nind>
    struct Numbers: NdField<specs::Number<T>, nind>{
        template<typename ...Args>
        Numbers(TableX *table, std::string description, Args... shape) :
            NdField<specs::Number<T>, nind>(table, {}, description, shape...) {}
    };

    template<typename T>
    using Number = Numbers<T, 0ul>;

    template<typename T, size_t nind, size_t nind_element>
    using NumberArrays = NdField<specs::NumberArray<T, nind_element>, nind>;

    template<size_t nind>
    using Bitsets = NdField<specs::Bitset, nind>;

    using Bitset = Bitsets<0ul>;

    template<size_t nind>
    using FermionOnvs = NdField<specs::FermionOnv, nind>;

    using FermionOnv = FermionOnvs<0ul>;

    template<size_t nind>
    using BosonOnvs = NdField<specs::BosonOnv, nind>;

    using BosonOnv = BosonOnvs<0ul>;

    template<size_t nind>
    using FermiBosOnvs = fb_onv::Field<nind>;

    using FermiBosOnv = FermiBosOnvs<0ul>;

    template<size_t nind, bool bosons>
    struct ConfigurationSelector {};

    template<size_t nind>
    struct ConfigurationSelector<nind, false> {
        typedef FermionOnvs<nind> type;
    };

    template<size_t nind>
    struct ConfigurationSelector<nind, true> {
        typedef FermiBosOnvs<nind> type;
    };

    template<size_t nind>
    using Onvs = typename ConfigurationSelector<nind, defs::bosons>::type;

    using Onv = Onvs<0ul>;
}

#endif //M7_FIELDS_H
