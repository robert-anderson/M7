//
// Created by anderson on 1/30/22.
//

#ifndef M7_XONVFIELD_H
#define M7_XONVFIELD_H

#include <M7_lib/basis/BasisData.h>

#include "CompositeField.h"
#include "Row.h"
#include "FrmOnvField.h"
#include "BosOnvField.h"
#include "FrmBosOnvField.h"

template<typename T>
struct XonvField : CompositeField<T, T> {
    using CompositeFieldBase::prefix;
    T m_ket, m_bra;

    XonvField(Row *row, HilbertSpace hs, std::string name) :
        CompositeField<T, T>(m_ket, m_bra),
        m_ket(row, hs, prefix("ket", name)),
        m_bra(row, hs, prefix("bra", name)) {}
};


struct FrmXonvField : XonvField<FrmOnvField> {
    FrmXonvField(Row *row, const HilbertSpace &hs, std::string name = "") :
            XonvField<FrmOnvField>(row, hs, name) {}

    FrmXonvField(Row *row, const FrmHilbertSpace &hs, std::string name = "") :
            FrmXonvField(row, HilbertSpace(hs, {0ul}), name) {}

    FrmXonvField &operator=(const std::pair<defs::inds, defs::inds> &inds) {
        m_ket = inds.first;
        m_bra = inds.second;
        return *this;
    }
};


struct BosXonvField : XonvField<BosOnvField> {
    BosXonvField(Row *row, const HilbertSpace &hs, std::string name = "") :
            XonvField<BosOnvField>(row, hs, name) {}

    BosXonvField(Row *row, const BosHilbertSpace &hs, std::string name = "") :
            BosXonvField(row, HilbertSpace({0ul, 0ul}, hs), name) {}

    BosXonvField &operator=(const std::pair<defs::inds, defs::inds> &inds) {
        m_ket = inds.first;
        m_bra = inds.second;
        return *this;
    }
};

struct FrmBosXonvField : XonvField<FrmBosOnvField> {
    FrmBosXonvField(Row *row, const HilbertSpace& hs, std::string name = "") :
        XonvField<FrmBosOnvField>(row, hs, name) {}

    FrmBosXonvField &
    operator=(const std::pair<std::pair<defs::inds, defs::inds>, std::pair<defs::inds, defs::inds>> &inds) {
        m_ket = inds.first;
        m_bra = inds.second;
        return *this;
    }
};

#endif //M7_XONVFIELD_H
