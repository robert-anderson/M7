//
// Created by Robert J. Anderson on 11/08/2021.
//

#ifndef M7_SPECTRALMOMENT_H
#define M7_SPECTRALMOMENT_H

#include <M7_lib/mae/MaeTable.h>
#include <M7_lib/conf/Conf.h>
#include <M7_lib/field/Fields.h>


//Communicator<MaeRow<wf_t>, MaeRow<wf_t>, true>
struct SpectralMoment {
    const OpSig m_exsig;
    const uint_t m_order;
    SpectralMoment(OpSig exsig, uint_t order): m_exsig(exsig), m_order(order){
        REQUIRE_EQ(order, 1ul, "Spectral moment quantities are currently only implemented for n=1")
        REQUIRE_TRUE(m_exsig.is_pure_frm(), "Spectral moment excitations must refer to purely fermionic perturbations");
    }
};

class SpecMoms {

public:

    SpecMoms(const conf::SpecMoms& /*opts*/) {}

    operator bool() const {
        return false;//!m_active_ranksigs.empty();
    }

    bool all_stores_empty() const {
//        for (auto& ranksig: m_active_ranksigs)
//            if (!m_rdms[ranksig]->m_store.is_freed())
//                return false;
            return true;
    }
};

#endif //M7_SPECTRALMOMENT_H
