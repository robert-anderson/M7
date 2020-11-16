//
// Created by rja on 04/06/2020.
//

#ifndef M7_EXCITATIONGENERATOR_H
#define M7_EXCITATIONGENERATOR_H

#include <src/core/hamiltonian/Hamiltonian.h>
#include <src/core/basis/Connections.h>
#include "PRNG.h"

class ExcitationGenerator {
protected:
    const Hamiltonian *m_h;
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

    virtual bool draw_single(const views::FermionOnv &src_fonv, views::Onv &dst_fonv,
                     const OccupiedOrbitals &occ, const VacantOrbitals &vac,
                     defs::prob_t &prob, defs::ham_t &helem, conn::AsFermionOnv &anticonn) = 0;

    virtual bool draw_double(const views::FermionOnv &src_fonv, views::FermionOnv &dst_fonv,
                     const OccupiedOrbitals &occ, defs::prob_t &prob, defs::ham_t &helem,
                             conn::AsFermionOnv &anticonn) = 0;

    const Hamiltonian *ham(){return m_h;}
};


#endif //M7_EXCITATIONGENERATOR_H
