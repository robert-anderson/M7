//
// Created by Robert J. Anderson on 10/08/2021.
//

#ifndef M7_FRMBOSONVCONNECTION_H
#define M7_FRMBOSONVCONNECTION_H

#include <M7_lib/field/FrmBosOnvField.h>

#include "FrmOnvConnection.h"
#include "BosOnvConnection.h"
#include "ComOps.h"

struct FrmBosOnvConnection {
    FrmOnvConnection m_frm;
    BosOnvConnection m_bos;

    explicit FrmBosOnvConnection(sys::Size size);

    explicit FrmBosOnvConnection(const FrmBosOnvField& mbf);

    void clear();

    void connect(const FrmBosOnvField& src, const FrmBosOnvField& dst);

    bool connect(const FrmBosOnvField &src, const FrmBosOnvField &dst, com_ops::FrmBos &com);

    void apply(const FrmBosOnvField& src, FrmBosOnvField& dst) const;

    size_t exsig() const;

    bool respects_occ_range(const FrmBosOnvField& src, size_t nboson_max) const;
};


#endif //M7_FRMBOSONVCONNECTION_H
