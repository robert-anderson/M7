//
// Created by rja on 10/08/2021.
//

#ifndef M7_FRMBOSONVCONNECTION_H
#define M7_FRMBOSONVCONNECTION_H

#include "FrmOnvConnection.h"
#include "BosOnvConnection.h"

struct FrmBosOnvConnection {
    FrmOnvConnection m_frm;
    BosOnvConnection m_bos;
    FrmBosOnvConnection(size_t nsite);

    void clear();

    void connect(const fields::FrmBosOnv& src, const fields::FrmBosOnv& dst);

    void apply(const fields::FrmBosOnv& src, fields::FrmBosOnv& dst) const;

    size_t exsig() const;
};


#endif //M7_FRMBOSONVCONNECTION_H
