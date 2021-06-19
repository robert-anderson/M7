//
// Created by jhalson on 28/04/2021.
//

#ifndef M7_SIGNPROBLEMFREETWF_H
#define M7_SIGNPROBLEMFREETWF_H

#include <src/core/enumerator/HamiltonianConnectionEnumerator.h>
#include <src/core/parallel/Reducible.h>
#include "src/core/dynamics/WalkerTable.h"
#include "src/core/hamiltonian/Hamiltonian.h"

class SpfTwfBase {
protected:
    std::vector<defs::ham_t> m_numerator;
    std::vector<defs::ham_t> m_denominator;
    size_t m_nsite;

public:
    std::vector<defs::ham_t> m_numerator_total;
    std::vector<defs::ham_t> m_denominator_total;

    SpfTwfBase(size_t npart, size_t nsite) :
            m_numerator(npart, 0.0), m_denominator(npart, 0.0),
            m_nsite(nsite),
            m_numerator_total(npart, 0.0), m_denominator_total(npart, 0.0) {
    };

    virtual void add(const Hamiltonian<0> &ham,
                     const fields::Numbers<defs::wf_t, defs::ndim_wf> &weight,
                     const fields::Onv<0> &onv) = 0;

    virtual void add(const Hamiltonian<1> &ham,
                     const fields::Numbers<defs::wf_t, defs::ndim_wf> &weight,
                     const fields::Onv<1> &onv) = 0;

    virtual void reduce() {
        mpi::all_sum(m_numerator.data(), m_numerator_total.data(), m_numerator.size());
        mpi::all_sum(m_denominator.data(), m_denominator_total.data(), m_denominator.size());
        m_numerator.assign(m_numerator.size(), 0.0);
        m_denominator.assign(m_denominator.size(), 0.0);
    };
};

#endif //M7_SIGNPROBLEMFREETWF_H
