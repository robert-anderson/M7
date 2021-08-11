//
// Created by rja on 11/08/2021.
//

#ifndef M7_FRMBOSONVFIELD_H
#define M7_FRMBOSONVFIELD_H

#include "FrmOnvField.h"
#include "BosOnvField.h"

struct FrmBosOnvField : MultiField<FrmOnvField, BosOnvField> {
    const std::string m_name;
    FrmOnvField &m_frm;
    BosOnvField &m_bos;

    FrmBosOnvField(Row *row, size_t nsite, std::string name = "");

    FrmBosOnvField(const FrmBosOnvField &other);

    FrmBosOnvField &operator=(const FrmBosOnvField &other);

    FrmBosOnvField &operator=(const std::pair<defs::inds, defs::inds> &inds);

    const size_t& nsite() const;
};


#endif //M7_FRMBOSONVFIELD_H
