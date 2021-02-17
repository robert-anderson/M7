//
// Created by rja on 09/02/2021.
//

#ifndef M7_FIELDSZ_H
#define M7_FIELDSZ_H

#include "NumberFieldZ.h"
#include "FermionOnvFieldZ.h"
#include "FlagFieldZ.h"
#include "NdMultiFieldZ.h"


/*
 * The field types can be summarized as follows:
 *
 *                items   elements
 * _________________________________
 * Number        |  1         1
 * Vector        |  1        many
 * Vectors       | many      many
 * FermionOnv    |  1        many
 * FermionOnvs   | many      many
 * Flag          |  1         1
 * Flags         |  1        many
 * FermiBosOnv   |  1        many
 * FermiBosOnvs  | many      many
 */

namespace fieldsz {

    template<typename T>
    struct Number : FieldZ<NumberFieldZ<T>> {
        Number(RowZ* row): FieldZ<NumberFieldZ<T>>(row, {}){}
    };

    template<typename T>
    struct Numbers : FieldZ<VectorFieldZ<T>> {
        Numbers(RowZ* row, size_t nelement): FieldZ<VectorFieldZ<T>>(row, {nelement}){}
    };

    template<typename T>
    struct NumberVectors : FieldZ<VectorsFieldZ<T>> {
        NumberVectors(RowZ* row, size_t nitem, size_t nelement): FieldZ<VectorsFieldZ<T>>(row, {nitem, nelement}){}
    };

    struct FermionOnv : FieldZ<FermionOnvFieldZ> {
        FermionOnv(RowZ* row, size_t nsite): FieldZ<FermionOnvFieldZ>(row, {nsite}){}
    };

    struct FermionOnvs : FieldZ<FermionOnvsFieldZ> {
        FermionOnvs(RowZ* row, size_t nitem, size_t nsite): FieldZ<FermionOnvsFieldZ>(row, {nitem, nsite}){}
    };

    struct Flag : FieldZ<FlagFieldZ> {
        Flag(RowZ* row): FieldZ<FlagFieldZ>(row, {}){}
    };

    struct Flags : FieldZ<FlagsFieldZ> {
        Flags(RowZ* row, size_t nflag): FieldZ<FlagsFieldZ>(row, {nflag}){}
    };

    using BosonOnvFieldZ = VectorFieldZ<uint8_t>;
    using BosonOnvsFieldZ = VectorsFieldZ<uint8_t>;

    struct FermiBosOnv : MultiFieldZ<FermionOnvFieldZ, BosonOnvFieldZ> {
        FermionOnvFieldZ& m_fonv;
        BosonOnvFieldZ& m_bonv;
        FermiBosOnv(RowZ* row, size_t nsite) :
                MultiFieldZ<FermionOnvFieldZ, BosonOnvFieldZ>(row, {nsite}, {nsite}),
                m_fonv(std::get<0>(m_subfields)), m_bonv(std::get<1>(m_subfields)){}
    };

    struct FermiBosOnvs : MultiFieldZ<FermionOnvsFieldZ, BosonOnvsFieldZ> {
        FermionOnvsFieldZ& m_fonv;
        BosonOnvsFieldZ& m_bonv;
        FermiBosOnvs(RowZ* row, size_t nitem, size_t nsite) :
                MultiFieldZ<FermionOnvsFieldZ, BosonOnvsFieldZ>(row, {nitem, nsite}, {nitem, nsite}),
                m_fonv(std::get<0>(m_subfields)), m_bonv(std::get<1>(m_subfields)){}
    };

    template<bool enable_bosons=defs::enable_bosons>
    using Onvs = typename std::conditional<enable_bosons, FermiBosOnvs, FermionOnvs>::type;

    template<bool enable_bosons=defs::enable_bosons>
    using Onv = typename std::conditional<enable_bosons, FermiBosOnv, FermionOnv>::type;

}


#endif //M7_FIELDSZ_H
