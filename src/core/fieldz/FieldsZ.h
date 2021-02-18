//
// Created by rja on 09/02/2021.
//

#ifndef M7_FIELDSZ_H
#define M7_FIELDSZ_H

#include "NumberFieldZ.h"
#include "FermionOnvFieldZ.h"
#include "FlagFieldZ.h"
#include "MultiFieldZ.h"


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

    struct BosonOnv : VectorFieldZ<uint8_t> {
        BosonOnv(RowZ *row, size_t nmode) : VectorFieldZ<uint8_t>(row, nmode) {
            ASSERT(m_nelement==nmode);
        }

        BosonOnv(const BosonOnv &other) : BosonOnv(other.m_row ? other.m_row->m_child : nullptr, other.m_nelement) {}

        BosonOnv &operator=(const defs::inds &inds) {
            MPI_ASSERT(inds.size() == m_nelement, "Vector is not the correct size");
            for (size_t i = 0ul; i < inds.size(); ++i) this->operator()(i) = inds[i];
            return *this;
        }
    };

    using BosonOnvs = VectorsFieldZ<uint8_t>;

    struct FermiBosOnv : MultiFieldZ<FermionOnv, BosonOnv> {
        FermionOnv &m_fonv;
        BosonOnv &m_bonv;

        FermiBosOnv(RowZ *row, size_t nsite) :
                MultiFieldZ<FermionOnv, BosonOnv>(row, {nullptr, nsite}, {nullptr, nsite}),
                m_fonv(get<0>()), m_bonv(get<1>()) {
        }

        FermiBosOnv(const FermiBosOnv& other): FermiBosOnv(other.m_fonv.row_of_copy(), other.m_fonv.m_nsite){}

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

    struct FermiBosOnvs : MultiFieldZ<FermionOnvs, BosonOnvs> {
        FermionOnvs &m_fonv;
        BosonOnvs &m_bonv;

        FermiBosOnvs(RowZ *row, size_t nitem, size_t nsite) :
                MultiFieldZ<FermionOnvs, BosonOnvs>(row, {nullptr, nitem, nsite}, {nullptr, nitem, nsite}),
                m_fonv(std::get<0>(m_subfields)), m_bonv(std::get<1>(m_subfields)) {}
    };

    template<bool enable_bosons = defs::enable_bosons>
    using Onv = typename std::conditional<enable_bosons, FermiBosOnv, FermionOnv>::type;

    template<bool enable_bosons = defs::enable_bosons>
    using Onvs = typename std::conditional<enable_bosons, FermiBosOnvs, FermionOnvs>::type;

}


#endif //M7_FIELDSZ_H