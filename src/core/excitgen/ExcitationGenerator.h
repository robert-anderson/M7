//
// Created by rja on 04/06/2020.
//

#ifndef M7_EXCITATIONGENERATOR_H
#define M7_EXCITATIONGENERATOR_H

#include <src/core/hamiltonian/Hamiltonian.h>
#include <src/core/connection/Connections.h>
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
        log::info("Excitation generator base initialized");
    }

    virtual bool draw(const fields::FrmOnv &src_onv, fields::FrmOnv &dst_onv,
                      const OccupiedOrbitals &occs, const VacantOrbitals &vacs,
                      defs::prob_t &prob, defs::ham_t &helem, conn::FrmOnv &conn) = 0;

    const Hamiltonian *ham(){return m_h;}

    virtual ~ExcitationGenerator(){}
};


#endif //M7_EXCITATIONGENERATOR_H
