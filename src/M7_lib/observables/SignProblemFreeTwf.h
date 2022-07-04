//
// Created by jhalson on 28/04/2021.
//

#ifndef M7_SIGNPROBLEMFREETWF_H
#define M7_SIGNPROBLEMFREETWF_H

#include "M7_lib/wavefunction/WalkerTable.h"
#include "M7_lib/foreach/ConnForeachGroup.h"

class SpfTwfBase {
protected:
    const Hamiltonian& m_ham;
    ConnForeachGroup m_conn_iters;
    v_t<ham_t> m_numerator;
    v_t<ham_t> m_denominator;

public:
    v_t<ham_t> m_numerator_total;
    v_t<ham_t> m_denominator_total;

    SpfTwfBase(const Hamiltonian &ham, uint_t npart);

    virtual void add(const field::Numbers<wf_t, c_ndim_wf> &weight,
                     const field::FrmOnv &onv) = 0;

    virtual void add(const field::Numbers<wf_t, c_ndim_wf> &weight,
                     const field::FrmBosOnv &onv) = 0;

    virtual void add(const field::Numbers<wf_t, c_ndim_wf> &weight,
                     const field::BosOnv &onv) = 0;

    virtual void reduce();
};

#endif //M7_SIGNPROBLEMFREETWF_H
