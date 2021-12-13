//
// Created by rja on 24/11/2021.
//

#ifndef M7_EXCITITERS_HUBBARDSINGLES_H
#define M7_EXCITITERS_HUBBARDSINGLES_H

#include "ExcitIters.h"

namespace excititers {

    struct HubbardSingles : public Frm {
        using Frm::foreach;
        const bool m_pbc;

        HubbardSingles(const Hamiltonian &ham);

        void foreach(const field::FrmOnv &src, conn::FrmOnv &conn, const fn_c_t<field::FrmOnv> &body) override;

    };
}


#endif //M7_EXCITITERS_HUBBARDSINGLES_H
