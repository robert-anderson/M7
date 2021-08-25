//
// Created by rja on 25/08/2021.
//

#ifndef M7_HUBBARD1DSINGLES_H
#define M7_HUBBARD1DSINGLES_H

#include "FrmConserve.h"

namespace excititers {
    struct Hubbard1dSingles : public FrmConserve {
        using FrmConserve::foreach;
        const bool m_pbc;

        Hubbard1dSingles(const Hamiltonian &ham);

        void foreach(const field::FrmOnv &src, conn::FrmOnv &conn, const fn_c_t <field::FrmOnv> &body) override;

    };
}


#endif //M7_HUBBARD1DSINGLES_H
