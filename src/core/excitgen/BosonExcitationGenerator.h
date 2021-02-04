//
// Created by rja on 12/11/2020.
//

#ifndef M7_BOSONEXCITATIONGENERATOR_H
#define M7_BOSONEXCITATIONGENERATOR_H

#include "src/core/hamiltonian/BosonCouplings.h"
#include "ExcitationGenerator.h"

class BosonExcitationGenerator : public ExcitationGenerator {

    defs::ham_t get_helement(const Hamiltonian<0>* ham,
            const size_t& p, const size_t& q, const size_t& imode){
        return 0.0;
    }
    defs::ham_t get_helement(const Hamiltonian<1>* ham,
            const size_t& p, const size_t& q, const size_t& imode){
        return ham->bc().get_element_1(p, q, imode);
    }

protected:
    size_t m_nboson_max;
public:

    BosonExcitationGenerator(const FermiBosHamiltonian *ham, PRNG& prng, size_t nboson_max):
        ExcitationGenerator(ham, prng), m_nboson_max(nboson_max){}

    bool draw_(const views::Onv<0> &src_onv, views::Onv<0> &dst_onv, const OccupiedOrbitals &occs,
              const VacantOrbitals &vacs, defs::prob_t &prob, defs::ham_t &helem,
              conn::Antisym<0> &anticonn) {
        return false;
    }

    bool draw_(const views::Onv<1> &src_onv, views::Onv<1> &dst_onv, const OccupiedOrbitals &occs,
              const VacantOrbitals &vacs, defs::prob_t &prob, defs::ham_t &helem,
              conn::Antisym<1> &anticonn) {
        if(m_nboson_max == 0) return false;

#ifndef NDEBUG
        auto nmode = src_onv.m_bonv.nmode();
        ASSERT(dst_onv.m_bonv.nmode() == nmode)
        ASSERT(nmode == src_onv.m_fonv.nsite() and nmode == dst_onv.m_fonv.nsite())
#endif

        auto imode_excit = occs[m_prng.draw_uint(occs.size())] % src_onv.m_fonv.nsite();
        int change;
        auto curr_occ = src_onv.m_bonv(imode_excit);

        prob = 1.0/occs.size();
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
        helem = get_helement(m_h, imode_excit, imode_excit, com);
        return true;
    }

    bool draw(const views::Onv<> &src_onv, views::Onv<> &dst_onv, const OccupiedOrbitals &occs,
               const VacantOrbitals &vacs, defs::prob_t &prob, defs::ham_t &helem,
               conn::Antisym<> &anticonn) override {
        return draw_(src_onv, dst_onv, occs, vacs, prob, helem, anticonn);
    }

};


#endif //M7_BOSONEXCITATIONGENERATOR_H
