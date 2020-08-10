//
// Created by rja on 04/06/2020.
//

#ifndef M7_EXCITATIONGENERATOR_H
#define M7_EXCITATIONGENERATOR_H


#include <src/core/hamiltonian/Hamiltonian.h>
#include "PRNG.h"

class ExcitationGenerator {
protected:
    const Hamiltonian *m_h;
    PRNG &m_prng;
    const size_t m_norb;
    const size_t m_nelec;
    const size_t m_norb_pair;
    const size_t m_nelec_pair;
    const bool m_spin_conserving;
public:
    ExcitationGenerator(const Hamiltonian *h, PRNG &prng) :
            m_h(h), m_prng(prng),
            m_norb(m_h->nsite() * 2),
            m_nelec(m_h->nelec()),
            m_norb_pair(integer_utils::nspair(m_norb)),
            m_nelec_pair(integer_utils::nspair(m_nelec)),
            m_spin_conserving(m_h->spin_conserving()) {
        std::cout << "Excitation generator base initialized" << std::endl;
    }

    virtual bool draw_single(const DeterminantElement &src_det, DeterminantElement &dst_det,
                     const OccupiedOrbitals &occ, const VacantOrbitals &vac,
                     defs::prob_t &prob, defs::ham_t &helem, AntisymConnection &anticonn) = 0;

    virtual bool draw_double(const DeterminantElement &src_det, DeterminantElement &dst_det,
                     const OccupiedOrbitals &occ, defs::prob_t &prob, defs::ham_t &helem,
                     AntisymConnection &anticonn) = 0;

    const Hamiltonian *ham(){return m_h;}
};


#endif //M7_EXCITATIONGENERATOR_H
