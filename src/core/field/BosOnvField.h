//
// Created by rja on 07/04/2021.
//

#ifndef M7_BOSONVFIELD_H
#define M7_BOSONVFIELD_H

#include "NumberField.h"
#include "src/core/basis/BasisData.h"

struct BosOnvField : NdNumberField<defs::bos_occ_t, 1> {
    typedef NdNumberField<defs::bos_occ_t, 1> base_t;
    /**
     * alias for the number of elements in the 1D numeric array
     */
    const size_t& m_nmode;

    BosOnvField(Row *row, size_t nmode, std::string name = "");

    BosOnvField(Row *row, BasisData bd, std::string name = "");

    BosOnvField(const BosOnvField &other) : base_t(other), m_nmode(m_nelement) {}

    BosOnvField &operator=(const BosOnvField &other) {
        base_t::operator=(other);
        return *this;
    }

    BosOnvField &operator=(const defs::inds &inds);

};


#endif //M7_BOSONVFIELD_H
