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
     * alias for the number of elements in the 1D numeric array
     */
    const BosHilbertSpace m_hs;
    /**
     * a refreshable cache of useful representations for excitation generation and enumeration
     */
    mutable decoded_mbf::BosOnv m_decoded;

    BosOnvField(Row *row, const BosHilbertSpace& hs, std::string name = "");

    BosOnvField(Row *row, const HilbertSpace& hs, std::string name = "");

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

};


#endif //M7_BOSONVFIELD_H
