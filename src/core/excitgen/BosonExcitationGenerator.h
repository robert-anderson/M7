//
// Created by rja on 12/11/2020.
//

#ifndef M7_BOSONEXCITATIONGENERATOR_H
#define M7_BOSONEXCITATIONGENERATOR_H

#include "src/core/hamiltonian/BosonCouplings.h"
#include "ExcitationGenerator.h"

#if 0
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

    bool draw_(const fields::Onv<0> &src_onv, fields::Onv<0> &dst_onv, const OccupiedOrbitals &occs,
               const VacantOrbitals &vacs, defs::prob_t &prob, defs::ham_t &helem,
               conn::Antisym<0> &anticonn) {
        return false;
    }

    bool draw_(const fields::Onv<1> &src_onv, fields::Onv<1> &dst_onv, const OccupiedOrbitals &occs,
               const VacantOrbitals &vacs, defs::prob_t &prob, defs::ham_t &helem,
               conn::Antisym<1> &anticonn) {
        if(m_nboson_max == 0) return false;

#ifndef NDEBUG
        auto nmode = src_onv.m_bos.nelement();
        ASSERT(nmode == src_onv.m_frm.m_nsite and nmode == dst_onv.m_frm.m_nsite)
#endif

        auto imode_excit = occs[m_prng.draw_uint(occs.size())] % src_onv.m_frm.m_nsite;
        int change;
        auto curr_occ = src_onv.m_bos[imode_excit];

        prob = 1.0/occs.size();

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

        auto com = src_onv.m_bos[imode_excit];
        if (change<0) com+=change;
        helem = get_helement(m_h, imode_excit, imode_excit, com);
        return true;
    }

    bool draw(const fields::Onv<> &src_onv, fields::Onv<> &dst_onv, const OccupiedOrbitals &occs,
              const VacantOrbitals &vacs, defs::prob_t &prob, defs::ham_t &helem,
              conn::Antisym<> &anticonn) override {
        return draw_(src_onv, dst_onv, occs, vacs, prob, helem, anticonn);
    }

};


#endif //M7_BOSONEXCITATIONGENERATOR_H
#endif //M7_BOSONEXCITATIONGENERATOR_H
