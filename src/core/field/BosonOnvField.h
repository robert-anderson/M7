//
// Created by rja on 07/04/2021.
//

#ifndef M7_BOSONONVFIELD_H
#define M7_BOSONONVFIELD_H

#include "NumberField.h"

struct BosonOnvField : NdNumberField<uint8_t, 1> {
    BosonOnvField(Row *row, size_t nmode) : NdNumberField<uint8_t, 1>(row, {nmode}) {}

    BosonOnvField(const BosonOnvField &other): BosonOnvField(other.m_row, other.m_format.extent(0)){}

    BosonOnvField &operator=(const BosonOnvField &other) {
        NdNumberField<uint8_t, 1>::operator=(other);
        return *this;
    }

    BosonOnvField &operator=(const defs::inds &inds) {
        MPI_ASSERT(inds.size() == nelement(), "Vector is not the correct size");
        for (size_t i = 0ul; i < inds.size(); ++i) this->operator[](i) = inds[i];
        return *this;
    }
};


#endif //M7_BOSONONVFIELD_H
