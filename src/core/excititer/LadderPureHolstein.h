//
// Created by rja on 24/11/2021.
//

#ifndef M7_LADDERPUREHOLSTEIN_H
#define M7_LADDERPUREHOLSTEIN_H

#include "ExcitIters.h"

namespace excititers {
    struct LadderPureHolstein : public Ladder {
        using Ladder::foreach;

        LadderPureHolstein(const Hamiltonian &ham, size_t exsig);

        void foreach(const field::FrmBosOnv &src, conn::FrmBosOnv &conn, const fn_c_t<field::FrmBosOnv> &body) override;
    };
}


#endif //M7_LADDERPUREHOLSTEIN_H
