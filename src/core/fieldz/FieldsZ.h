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
    using Number = NumberFieldZ<T>;

    template<typename T>
    using Vector = VectorFieldZ<T>;

    template<typename T>
    using Vectors = VectorsFieldZ<T>;

    using FermionOnv = FermionOnvFieldZ;
    using FermionOnvs = FermionOnvsFieldZ;

    using Flag = FlagFieldZ;
    using Flags = FlagsFieldZ;

    using BosonOnv = VectorFieldZ<uint8_t>;
    using BosonOnvs = VectorsFieldZ<uint8_t>;

    struct FermiBosOnv : MultiFieldZ<FermionOnv, BosonOnv> {
        FermionOnv& m_fonv;
        BosonOnv& m_bonv;
        FermiBosOnv(RowZ* row, size_t nsite) :
                MultiFieldZ<FermionOnv, BosonOnv>(row, {nullptr, nsite}, {nullptr, nsite}),
                m_fonv(std::get<0>(m_subfields)), m_bonv(std::get<1>(m_subfields)){}
    };

    struct FermiBosOnvs : MultiFieldZ<FermionOnvs, BosonOnvs> {
        FermionOnvs& m_fonv;
        BosonOnvs& m_bonv;
        FermiBosOnvs(RowZ* row, size_t nitem, size_t nsite) :
                MultiFieldZ<FermionOnvs, BosonOnvs>(row, {nullptr, nitem, nsite}, {nullptr, nitem, nsite}),
                m_fonv(std::get<0>(m_subfields)), m_bonv(std::get<1>(m_subfields)){}
    };

    template<bool enable_bosons=defs::enable_bosons>
    using Onv = typename std::conditional<enable_bosons, FermiBosOnv, FermionOnv>::type;

    template<bool enable_bosons=defs::enable_bosons>
    using Onvs = typename std::conditional<enable_bosons, FermiBosOnvs, FermionOnvs>::type;

}


#endif //M7_FIELDSZ_H
