//
// Created by rja on 25/08/2021.
//

#ifndef M7_FRMCONSERVE_H
#define M7_FRMCONSERVE_H

#include <src/core/util/Foreach.h>
#include "ExcitIter.h"

namespace excititers {

    struct FrmConserve : public Frm {
        using Frm::foreach;
    protected:
        foreach::rtnd::Ordered<> m_cre_loop;
        foreach::rtnd::Ordered<> m_ann_loop;

        fn_c_t <field::FrmOnv> convert(conn::FrmBosOnv &work_conn, const fn_c_t <field::FrmBosOnv> &fn);

    public:
        FrmConserve(const Hamiltonian &ham, size_t exsig);

        void foreach(const field::FrmOnv &src, conn::FrmOnv &conn, const fn_c_t <field::FrmOnv> &body) override;

        void foreach(const FrmBosOnv &src, conn::FrmBosOnv &conn, const fn_c_t <FrmBosOnv> &body) override;

        void foreach(const BosOnv &src, conn::BosOnv &conn, const fn_c_t <BosOnv> &body) override {}
    };

}

#endif //M7_FRMCONSERVE_H
