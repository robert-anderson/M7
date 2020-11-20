//
// Created by rja on 04/06/2020.
//

#ifndef M7_EXCITATIONGENERATOR_H
#define M7_EXCITATIONGENERATOR_H

#include <src/core/hamiltonian/Hamiltonian.h>
#include <src/core/basis/Connections.h>
#include "src/core/sample/PRNG.h"

class ExcitationGenerator {
protected:
    const Hamiltonian *m_h;
    PRNG &m_prng;
    const size_t m_nintind;
    const size_t m_nelec;
    const size_t m_norb_pair;
    const size_t m_nelec_pair;
public:
    ExcitationGenerator(const Hamiltonian *h, PRNG &prng) :
            m_h(h), m_prng(prng),
            m_nintind(m_h->nsite() * 2),
            m_nelec(m_h->nelec()),
            m_norb_pair(integer_utils::nspair(m_nintind)),
            m_nelec_pair(integer_utils::nspair(m_nelec)) {
        std::cout << "Excitation generator base initialized" << std::endl;
    }

    virtual bool draw(const views::FermionOnv &src_fonv, views::FermionOnv &dst_fonv,
                            const OccupiedOrbitals &occ, const VacantOrbitals &vac,
                            defs::prob_t &prob, defs::ham_t &helem, conn::AsFermionOnv &anticonn) = 0;

    virtual bool draw(const views::FermiBosOnv &src_fonv, views::FermiBosOnv &dst_fonv,
                            const OccupiedOrbitals &occ, const VacantOrbitals &vac,
                            defs::prob_t &prob, defs::ham_t &helem, conn::AsFermiBosOnv &anticonn) = 0;

    const Hamiltonian *ham(){return m_h;}
};


#endif //M7_EXCITATIONGENERATOR_H
