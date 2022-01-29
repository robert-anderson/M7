//
// Created by rja on 29/01/2022.
//

#ifndef M7_FRMXONVFIELD_H
#define M7_FRMXONVFIELD_H

#include "CompositeField.h"
#include "FrmOnvField.h"

struct FrmXonvField : CompositeField<FrmOnvField, FrmOnvField> {
    FrmOnvField m_ket, m_bra;

    FrmXonvField(Row* row, size_t nsite, std::string name="");

    FrmXonvField(Row* row, BasisDims bd, std::string name="");

    FrmXonvField(const FrmXonvField& other);

    FrmXonvField& operator=(const FrmXonvField& other) {
        m_ket = other.m_ket;
        m_bra = other.m_bra;
        return *this;
    }

    FrmXonvField& operator=(const std::pair<defs::inds, defs::inds> &inds);
};

#endif //M7_FRMXONVFIELD_H
