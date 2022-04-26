//
// Created by rja on 07/04/2021.
//

#ifndef M7_BOSONVFIELD_H
#define M7_BOSONVFIELD_H

#include <M7_lib/basis/BasisData.h>
#include <M7_lib/caches/DecodedBosOnv.h>

#include "NumberField.h"

struct BosOnvField : NdNumberField<defs::bos_occ_t, 1> {
    typedef NdNumberField<defs::bos_occ_t, 1> base_t;
    /**
     * single-particle basis data is all the system-related data that this class needs to retain
     */
    const sys::bos::Basis m_basis;
    /**
     * a refreshable cache of useful representations for excitation generation and enumeration
     */
    mutable decoded_mbf::BosOnv m_decoded;

    BosOnvField(Row *row, const sys::bos::Basis& basis, std::string name = "");
    /*
     * all ONVs implement the following ctor
     */
    BosOnvField(Row* row, const sys::Basis& basis, std::string name="");
    /*
     * this particular MBF only needs the basis, but future MBF types might need the full sector information, and so
     * a common interface is realised by implementing a ctor of the following form in all MBFs
     */
    BosOnvField(Row *row, const sys::Sector& hs, std::string name = "");

    BosOnvField(const BosOnvField &other);

    BosOnvField &operator=(const BosOnvField &other) {
        base_t::operator=(other);
        return *this;
    }

    BosOnvField &operator=(const defs::inds &inds);

    /**
     * set based on the (repeated) boson operator indices
     * @param iops
     *  operator indices (one for each boson in the state)
     */
    void set_ops(const defs::inds &iops);

    size_t nboson() const;

};


#endif //M7_BOSONVFIELD_H
