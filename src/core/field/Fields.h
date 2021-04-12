//
// Created by rja on 09/02/2021.
//

#ifndef M7_FIELDS_H
#define M7_FIELDS_H

#include "NumberField.h"
#include "MultiField.h"
#include "FermionOnvField.h"
#include "BosonOnvField.h"

namespace fields {
    template<typename T, size_t nind>
    using Numbers = NdNumberField<T, nind>;

    template<typename T>
    using Number = NumberField<T>;

    template<size_t nind>
    using Flags = BitsetField<uint8_t, nind>;
    using Flag = BitField<uint8_t>;

    using FermionOnv = FermionOnvField;

    using BosonOnv = BosonOnvField;

    struct FermiBosOnv : MultiField<FermionOnv, BosonOnv> {
        const std::string m_name;
        FermionOnv &m_fonv;
        BosonOnv &m_bonv;

        FermiBosOnv(Row *row, size_t nsite, std::string name="") :
                MultiField<FermionOnv, BosonOnv>(row,
                                                 {nullptr, nsite, name.empty() ? "" : name+" (fermion)"},
                                                 {nullptr, nsite, name.empty() ? "" : name+" (boson)"}),
                m_name(name), m_fonv(get<0>()), m_bonv(get<1>()) {
        }

        FermiBosOnv(const FermiBosOnv& other): FermiBosOnv(other.m_fonv.row_of_copy(), other.m_fonv.m_nsite, other.m_name){}

        FermiBosOnv &operator=(const FermiBosOnv& other) {
            m_fonv = other.m_fonv;
            m_bonv = other.m_bonv;
            return *this;
        }

        FermiBosOnv &operator=(const std::pair<defs::inds, defs::inds>& inds) {
            m_fonv = inds.first;
            m_bonv = inds.second;
            return *this;
        }
    };

    template<bool enable_bosons = defs::enable_bosons>
    using Onv = typename std::conditional<enable_bosons, FermiBosOnv, FermionOnv>::type;
}

#endif //M7_FIELDS_H