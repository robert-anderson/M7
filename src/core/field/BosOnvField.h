//
// Created by rja on 07/04/2021.
//

#ifndef M7_BOSONVFIELD_H
#define M7_BOSONVFIELD_H

#include "NumberField.h"
#include "src/core/basis/BasisDims.h"

struct BosOnvField : NdNumberField<defs::bos_occ_t, 1> {
    BosOnvField(Row *row, BasisDims bd, std::string name="") :
    NdNumberField<defs::bos_occ_t, 1>(row, {{bd.m_nmode}, {"boson mode occupations"}}, name) {
        bd.require_pure_bos();
    }

    BosOnvField(const BosOnvField &other): BosOnvField(other.m_row, {0ul, other.m_format.m_shape[0]}, other.m_name){}

    BosOnvField &operator=(const BosOnvField &other) {
        NdNumberField<defs::bos_occ_t, 1>::operator=(other);
        return *this;
    }

    BosOnvField &operator=(const defs::inds &inds) {
        DEBUG_ASSERT_EQ(inds.size(), nelement(), "Vector is not the correct size");
        for (size_t i = 0ul; i < inds.size(); ++i) this->operator[](i) = inds[i];
        return *this;
    }
};


#endif //M7_BOSONVFIELD_H
