//
// Created by jhalson on 28/04/2021.
//

#ifndef M7_SIGNPROBLEMFREETWF_H
#define M7_SIGNPROBLEMFREETWF_H

#include "src/core/wavefunction/WalkerTable.h"
#include "src/core/hamiltonian/ForeachConnection.h"

class SpfTwfBase {
protected:
    const Hamiltonian& m_ham;
    foreach_conn::vector_t m_foreach_conns;
    std::vector<defs::ham_t> m_numerator;
    std::vector<defs::ham_t> m_denominator;
    size_t m_nsite;

public:
    std::vector<defs::ham_t> m_numerator_total;
    std::vector<defs::ham_t> m_denominator_total;

    SpfTwfBase(const Hamiltonian &ham, size_t npart, size_t nsite) :
            m_ham(ham),  m_foreach_conns(foreach_conn::make_all(ham)),
            m_numerator(npart, 0.0), m_denominator(npart, 0.0),
            m_nsite(nsite), m_numerator_total(npart, 0.0), m_denominator_total(npart, 0.0) {
    }

    virtual void add(const field::Numbers<defs::wf_t, defs::ndim_wf> &weight,
                     const field::FrmOnv &onv) = 0;

    virtual void add(const field::Numbers<defs::wf_t, defs::ndim_wf> &weight,
                     const field::FrmBosOnv &onv) = 0;

    virtual void reduce() {
        mpi::all_sum(m_numerator.data(), m_numerator_total.data(), m_numerator.size());
        mpi::all_sum(m_denominator.data(), m_denominator_total.data(), m_denominator.size());
        m_numerator.assign(m_numerator.size(), 0.0);
        m_denominator.assign(m_denominator.size(), 0.0);
    };
};

#endif //M7_SIGNPROBLEMFREETWF_H
