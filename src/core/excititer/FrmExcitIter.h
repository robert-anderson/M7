//
// Created by rja on 25/08/2021.
//

#ifndef M7_FRMEXCITITER_H
#define M7_FRMEXCITITER_H

#include <src/core/util/Foreach.h>
#include "ExcitIter.h"

struct FrmExcitIter : ExcitIter {
    using ExcitIter::foreach;
protected:
    foreach::rtnd::Ordered<> m_cre_loop;
    foreach::rtnd::Ordered<> m_ann_loop;

    fn_c_t<field::FrmOnv> convert(conn::FrmBosOnv &work_conn, const fn_c_t<field::FrmBosOnv>& fn);

public:
    FrmExcitIter(size_t exsig, const Hamiltonian& ham);

    void foreach(const field::FrmOnv& src, conn::FrmOnv& work_conn, const fn_c_t<field::FrmOnv> &body) override;

    void foreach(const FrmBosOnv &src, conn::FrmBosOnv &work_conn, const fn_c_t<FrmBosOnv> &body) override;

    void foreach(const BosOnv &src, conn::BosOnv &work_conn, const fn_c_t<BosOnv> &body) override {}
};

#endif //M7_FRMEXCITITER_H
