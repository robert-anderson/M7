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
    struct NumberArrays : NdFieldZ<nind, NumberFieldZ<T, nind_element>> {
        NumberArrays(RowZ *row, std::array<size_t, nind> shape, std::array<size_t, nind_element> element_shape) :
                NdFieldZ<nind, NumberFieldZ<T, nind_element>>(row, shape, element_shape) {}
    };

    template<typename T, size_t nind_element>
    struct NumberArray : NdFieldZ<0ul, NumberFieldZ<T, nind_element>> {
        NumberArray(RowZ *row, std::array<size_t, nind_element> element_shape) :
                NdFieldZ<0ul, NumberFieldZ<T, nind_element>>(row, {}, element_shape) {}
    };

    template<typename T, size_t nind>
    struct Numbers : NdFieldZ<nind, NumberFieldZ<T, 0ul>> {
        Numbers(RowZ *row, std::array<size_t, nind> shape) :
                NdFieldZ<nind, NumberFieldZ<T, 0ul>>(row, shape, {}) {}
    };

    template<typename T>
    struct Number : NdFieldZ<0ul, NumberFieldZ<T, 0ul>> {
        Number(RowZ *row) : NdFieldZ<0ul, NumberFieldZ<T, 0ul>>(row, {}, {}) {}
    };

    template<size_t nind>
    struct BosonOnvs : NdFieldZ<nind, BosonOnvFieldZ> {
        BosonOnvs(RowZ *row, std::array<size_t, nind> shape, size_t nmode) :
                NdFieldZ<nind, BosonOnvFieldZ>(row, shape, {nmode}) {}
    };

    struct BosonOnv : BosonOnvs<0ul> {
        BosonOnv(RowZ *row, size_t nmode) : BosonOnvs<0ul>(row, {}, nmode) {}
    };

    template<size_t nind_element>
    struct Flags : NdFieldZ<0ul, FlagFieldZ<nind_element>> {
        Flags(RowZ *row, std::array<size_t, nind_element> element_shape) :
                NdFieldZ<0ul, FlagFieldZ<nind_element>>(row, {}, element_shape) {}
    };

    struct Flag : Flags<0ul> {
        Flag(RowZ *row) : Flags<0ul>(row, {}) {}
    };

    template<size_t nind>
    struct FermionOnvs : NdFieldZ<nind, FermionOnvFieldZ> {
        FermionOnvs(RowZ *row, std::array<size_t, nind> shape, size_t nsite) :
                NdFieldZ<nind, FermionOnvFieldZ>(row, shape, {nsite}) {}
    };

    struct FermionOnv : FermionOnvs<0ul> {
        FermionOnv(RowZ *row, size_t nsite) : FermionOnvs<0ul>(row, {}, nsite) {}
    };

    template<size_t nind>
    struct FermiBosOnvs : NdMultiFieldZ<nind, FermionOnvFieldZ, BosonOnvFieldZ> {
        FermiBosOnvs(RowZ *row, std::array<size_t, nind> shape, size_t nsite) :
                NdMultiFieldZ<nind, FermionOnvFieldZ, BosonOnvFieldZ>(row, shape, {nsite}, {nsite}) {}

        using NdMultiFieldZ<nind, FermionOnvFieldZ, BosonOnvFieldZ>::get;

        const FermionOnvFieldZ &get_fonv() const {
            return this->template get<0>();
        }
        const FermionOnvFieldZ &get_bonv() const {
            return this->template get<1>();
        }
    };

    struct FermiBosOnv : FermiBosOnvs<0ul> {
        FermiBosOnv(RowZ *row, size_t nsite) : FermiBosOnvs<0ul>(row, {}, nsite) {}
    };

    template<size_t nind, bool enable_bosons = defs::enable_bosons>
    using Onvs = typename std::conditional<enable_bosons, FermiBosOnvs<nind>, FermionOnvs<nind>>::type;

    template<bool enable_bosons = defs::enable_bosons>
    using Onv = typename std::conditional<enable_bosons, FermiBosOnv, FermionOnv>::type;

}


#endif //M7_FIELDSZ_H
