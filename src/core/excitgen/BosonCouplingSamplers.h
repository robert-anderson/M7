//
// Created by rja on 12/11/2020.
//

#ifndef M7_BOSONCOUPLINGSAMPLERS_H
#define M7_BOSONCOUPLINGSAMPLERS_H

#include <src/core/basis/DecodedDeterminant.h>
#include <src/core/basis/Connections.h>
#include "src/core/field/Views.h"
#include "src/core/sample/PRNG.h"

struct BosonCouplingSamplers {

    size_t m_nboson_max;
    PRNG m_prng;
    BosonCouplingSamplers(size_t nboson_max, PRNG prng):m_nboson_max(nboson_max), m_prng(prng){}

    bool draw_single(const views::FermiBosOnv &src_onv, views::FermiBosOnv &dst_onv,
                     const OccupiedOrbitals &occ, defs::prob_t &prob, defs::ham_t &helem, conn::AsFermiBosOnv &anticonn) {
        /*
         * TODO: James
         *  make sure to either prevent the draw of a boson onv which exceeds the cutoff, or at
         *  least return false in such cases
         */
        auto nmode = src_onv.m_bonv.nmode();
        ASSERT(dst_onv.nmode() == nmode);
        auto excit_ind = m_prng.draw_uint(nmode);
        int change;
        auto curr_occ = src_onv.m_bonv(excit_ind);
        if(curr_occ == m_nboson_max){
            change = -1;
        }
        else{
            change = m_prng.draw_uint(2) == 0 ? -1 : 1;
        }
        dst_onv.m_bonv(excit_ind) += change;

        auto nmax = 0ul;
        for(auto imode = 0ul; imode < nmode; ++imode){
            nmax += src_onv.m_bonv(imode) == m_nboson_max;
        }
        prob = 1.0/(2*nmode - nmax);

        // TODO need to do fermion part with existing excitgen.



        return false;
    }

};


#endif //M7_BOSONCOUPLINGSAMPLERS_H
