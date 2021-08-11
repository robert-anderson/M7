//
// Created by rja on 11/08/2021.
//

#ifndef M7_SPECTRALMOMENT_H
#define M7_SPECTRALMOMENT_H

#include "src/core/field/Fields.h"


Communicator<MaeRow<defs::wf_t>, MaeRow<defs::wf_t>, true>, Archivable

struct SpectralMoment {
    const size_t m_exsig, m_order;
    SpectralMoment(size_t exsig, size_t order): m_exsig(exsig), m_order(order){
        REQUIRE_EQ(order, 1ul, "Spectral moment quantities are currently only implemented for n=1")
        REQUIRE_TRUE(conn_utils::pure_frm(exsig),
                     "Spectral moment excitations must refer to purely fermionic perturbations");
    }
};

#endif //M7_SPECTRALMOMENT_H
