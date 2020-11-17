//
// Created by rja on 12/11/2020.
//

#ifndef M7_BOSONCOUPLINGSAMPLERS_H
#define M7_BOSONCOUPLINGSAMPLERS_H

#include <src/core/basis/DecodedDeterminant.h>
#include <src/core/basis/Connections.h>
#include "src/core/field/Views.h"
#include "src/core/sample/PRNG.h"
#include "src/core/hamiltonian/BosonCouplings.h"

struct BosonCouplingSamplers {

    size_t m_nboson_max;
    PRNG& m_prng;
    const BosonCouplings& m_bc;
    BosonCouplingSamplers(const BosonCouplings& bc, size_t nboson_max, PRNG& prng):
        m_nboson_max(nboson_max), m_prng(prng), m_bc(bc) {}

    bool draw_single(const views::FermiBosOnv &src_onv, views::FermiBosOnv &dst_onv,
                     const OccupiedOrbitals &occ, defs::prob_t &prob, defs::ham_t &helem, conn::AsFermiBosOnv &anticonn) {

        if(m_nboson_max == 0) return false;

        auto nmode = src_onv.m_bonv.nmode();
        ASSERT(dst_onv.m_bonv.nmode() == nmode)
        ASSERT(nmode == src_onv.m_fonv.nsite() and nmode == dst_onv.m_fonv.nsite())

        auto imode_excit = occ.m_inds[m_prng.draw_uint(occ.m_nind)] % src_onv.m_fonv.nsite();
        int change;
        auto curr_occ = src_onv.m_bonv(imode_excit);

        prob = 1.0/occ.m_nind;
        if(src_onv.m_fonv.get(0, imode_excit) and src_onv.m_fonv.get(1, imode_excit)){
            prob *= 2.0;
        }

        if(curr_occ == m_nboson_max){
            change = -1;
        }
        else if(curr_occ == 0){
            change = 1;
        }
        else{
            change = m_prng.draw_uint(2) == 0 ? -1 : 1; // TODO factor this out
            prob *= 0.5;
        }

        anticonn.m_bonvconn.add(imode_excit, change);
        anticonn.apply(src_onv, dst_onv);
        helem = m_bc.get_element_1(anticonn);
        return true;
    }

};


#endif //M7_BOSONCOUPLINGSAMPLERS_H
