//
// Created by rja on 04/06/2020.
//

#ifndef M7_EXCITATIONGENERATOR_H
#define M7_EXCITATIONGENERATOR_H

#if 0

#include <src/core/hamiltonian/FermionHamiltonian.h>
#include "PRNG.h"

class ExcitationGenerator {
protected:
    const FermionHamiltonian *m_h;
    PRNG &m_prng;
    const size_t m_nintind;
    const size_t m_nelec;
    const size_t m_norb_pair;
    const size_t m_nelec_pair;
    const bool m_spin_conserving_1e, m_spin_conserving_2e;
public:
    ExcitationGenerator(const FermionHamiltonian *h, PRNG &prng) :
            m_h(h), m_prng(prng),
            m_nintind(m_h->nsite() * 2),
            m_nelec(m_h->nelec()),
            m_norb_pair(integer_utils::nspair(m_nintind)),
            m_nelec_pair(integer_utils::nspair(m_nelec)),
            m_spin_conserving_1e(m_h->spin_conserving_1e()),
            m_spin_conserving_2e(m_h->spin_conserving_2e()) {
        std::cout << "Excitation generator base initialized" << std::endl;
    }

    virtual bool draw_single(const DeterminantElement &src_det, DeterminantElement &dst_det,
                     const OccupiedOrbitals &occ, const VacantOrbitals &vac,
                     defs::prob_t &prob, defs::ham_t &helem, AntisymFermionOnvConnection &anticonn) = 0;

    virtual bool draw_double(const DeterminantElement &src_det, DeterminantElement &dst_det,
                     const OccupiedOrbitals &occ, defs::prob_t &prob, defs::ham_t &helem,
                     AntisymFermionOnvConnection &anticonn) = 0;

    const FermionHamiltonian *ham(){return m_h;}
};


#endif //M7_EXCITATIONGENERATOR_H
#endif //M7_EXCITATIONGENERATOR_H
