//
// Created by rja on 12/11/2020.
//

#ifndef M7_BOSONEXCITATIONGENERATOR_H
#define M7_BOSONEXCITATIONGENERATOR_H

#include "src/core/hamiltonian/BosonCouplings.h"
#include "ExcitationGenerator.h"

class BosonExcitationGenerator : public ExcitationGenerator {

protected:
    size_t m_nboson_max;
public:

    BosonExcitationGenerator(const FermiBosHamiltonian *ham, PRNG& prng, size_t nboson_max):
        ExcitationGenerator(ham, prng), m_nboson_max(nboson_max){}

    bool draw(const views::FermionOnv &src_fonv, views::FermionOnv &dst_fonv, const OccupiedOrbitals &occ,
              const VacantOrbitals &vac, defs::prob_t &prob, defs::ham_t &helem,
              conn::AsFermionOnv &anticonn) override {
        return false;
    }

    bool draw(const views::FermiBosOnv &src_onv, views::FermiBosOnv &dst_onv, const OccupiedOrbitals &occ,
              const VacantOrbitals &vac, defs::prob_t &prob, defs::ham_t &helem,
              conn::AsFermiBosOnv &anticonn) override {
        if(m_nboson_max == 0) return false;

        auto nmode = src_onv.m_bonv.nmode();
        ASSERT(dst_onv.m_bonv.nmode() == nmode)
        ASSERT(nmode == src_onv.m_fonv.nsite() and nmode == dst_onv.m_fonv.nsite())

        auto imode_excit = occ.m_inds[m_prng.draw_uint(occ.m_nind)] % src_onv.m_fonv.nsite();
        int change;
        auto curr_occ = src_onv.m_bonv(imode_excit);

        prob = 1.0/occ.m_nind;
        // there are two ways to generate such connections, and they should be twice as likely
//        if(src_onv.m_fonv.get(0, imode_excit) and src_onv.m_fonv.get(1, imode_excit)){
//            prob *= 1.0;
//        }

        if(curr_occ == m_nboson_max){
            change = -1;
        }
        else if(curr_occ == 0){
            change = 1;
        }
        else{
            change = m_prng.draw_uint(2) == 0 ? -1 : 1; // TODO use one PRNG for both by drawing from [0, 2*occ.m_nind)
            prob *= 0.5;
        }

        anticonn.zero();
        anticonn.m_bonvconn.add(imode_excit, change);
        anticonn.apply(src_onv, dst_onv);

        auto com = src_onv.m_bonv(imode_excit);
        if (change<0) com+=change;
        helem = m_h->bc().get_element_1(imode_excit, imode_excit, com);
        return true;
    }

};


#endif //M7_BOSONEXCITATIONGENERATOR_H
