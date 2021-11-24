//
// Created by rja on 24/11/2021.
//

#ifndef M7_LADDERPURE_H
#define M7_LADDERPURE_H

#include "ExcitIters.h"

namespace excititers {
    struct LadderPure : public Ladder {
        using Ladder::foreach;

        LadderPure(const Hamiltonian &ham, size_t exsig);

        void foreach(const field::FrmBosOnv &src, conn::FrmBosOnv &conn, const fn_c_t<field::FrmBosOnv> &body) override;
    };
}

#endif //M7_LADDERPURE_H
