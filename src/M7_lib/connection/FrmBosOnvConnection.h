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


    void connect(const FrmOnvField& src, const FrmOnvField& dst);

    bool connect(const FrmOnvField &src, const FrmOnvField &dst, com_ops::Frm &com);

    void connect(const BosOnvField& src, const BosOnvField& dst);

    bool connect(const BosOnvField &src, const BosOnvField &dst, com_ops::Bos &com);

    void connect(const FrmBosOnvField& src, const FrmBosOnvField& dst);

    bool connect(const FrmBosOnvField &src, const FrmBosOnvField &dst, com_ops::FrmBos &com);

    void apply(const FrmBosOnvField& src, FrmBosOnvField& dst) const;

    OpSig exsig() const;

    bool respects_occ_range(const FrmBosOnvField& src, uint_t nboson_max) const;

    bool phase(const FrmBosOnvField &src) const {
        return m_frm.phase(src.m_frm);
    }
};

static std::ostream &operator<<(std::ostream &os, const FrmBosOnvConnection &conn) {
    os << conn.m_frm << conn.m_bos;
    return os;
}

#endif //M7_FRMBOSONVCONNECTION_H
