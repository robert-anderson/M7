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

    XonvField(Row *row, BasisData bd, std::string name) :
        CompositeField<T, T>(m_ket, m_bra),
        m_ket(row, bd, prefix("ket", name)),
        m_bra(row, bd, prefix("bra", name)) {}
};


struct FrmXonvField : XonvField<FrmOnvField> {
    FrmXonvField(Row *row, const BasisData &bd, std::string name = "") :
            XonvField<FrmOnvField>(row, bd, name) {}

    FrmXonvField(Row *row, const FrmBasisData &bd, std::string name = "") :
            FrmXonvField(row, {bd, {}}, name) {}

    FrmXonvField &operator=(const std::pair<defs::inds, defs::inds> &inds) {
        m_ket = inds.first;
        m_bra = inds.second;
        return *this;
    }
};


struct BosXonvField : XonvField<BosOnvField> {
    BosXonvField(Row *row, const BasisData &bd, std::string name = "") :
            XonvField<BosOnvField>(row, bd, name) {}

    BosXonvField(Row *row, const BosBasisData &bd, std::string name = "") :
            BosXonvField(row, {{}, bd}, name) {}

    BosXonvField &operator=(const std::pair<defs::inds, defs::inds> &inds) {
        m_ket = inds.first;
        m_bra = inds.second;
        return *this;
    }
};

struct FrmBosXonvField : XonvField<FrmBosOnvField> {
    FrmBosXonvField(Row *row, BasisData bd, std::string name = "") : XonvField<FrmBosOnvField>(row, bd, name) {}

    FrmBosXonvField &
    operator=(const std::pair<std::pair<defs::inds, defs::inds>, std::pair<defs::inds, defs::inds>> &inds) {
        m_ket = inds.first;
        m_bra = inds.second;
        return *this;
    }
};

#endif //M7_XONVFIELD_H
