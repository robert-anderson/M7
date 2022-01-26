//
// Created by anderson on 1/25/22.
//

#ifndef M7_FRMBOSXONVFIELD_H
#define M7_FRMBOSXONVFIELD_H

#include "FrmBosOnvField.h"

#if 0
struct FrmBosXonvField : MultiField<FrmBosOnvField, FrmBosOnvField> {
    const std::string m_name;
    FrmBosOnvField &m_ket;
    FrmBosOnvField &m_bra;

    //FrmBosOnvField(Row *row, BasisDims bd, std::string name = "");

    FrmBosXonvField(Row *row, BasisDims bd, std::string name = ""):
            MultiField<FrmBosOnvField, FrmBosOnvField>(row,
                 {nullptr, bd, name.empty() ? "" : name + " (ket)"},
                 {nullptr, bd, name.empty() ? "" : name + " (bra)"}),
            m_name(name), m_ket(get<1>()), m_bra(get<1>()) {
    }
};

#endif //M7_FRMBOSXONVFIELD_H
#endif //M7_FRMBOSXONVFIELD_H
