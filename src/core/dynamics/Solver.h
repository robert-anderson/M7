//
// Created by rja on 10/11/2020.
//

#ifndef M7_SOLVER_H
#define M7_SOLVER_H

#include <src/core/hamiltonian/Hamiltonian.h>
#include "Reference.h"
#include "src/core/table/CommunicatingPair.h"
#include "src/core/field/Fields.h"
#include "Propagator.h"

class Solver {

    Propagator &m_prop;
    const Options &m_opts;
    Wavefunction &m_wf;
    Reference m_reference;

public:
    // std::string name, size_t nbucket, fields::Onv::params_t onv_params, size_t nroot, size_t nreplica
    Solver(Propagator &prop, Wavefunction &wf, views::Onv ref_onv) :
            m_prop(prop),
            m_opts(prop.m_opts),
            m_wf(wf),
            m_reference(m_wf, m_prop.m_ham, ref_onv, m_opts)
            {}

    void execute(size_t niter=1){
        for(size_t i=0ul; i<niter; ++i){
            propagate();
            m_wf.m_spawn.communicate();
            annihilate();
        }
    }
    void preloop() {
        //if (m_prop->semi_stochastic()) m_detsub->gather_and_project();
    }

    void propagate();

    void annihilate_row(const size_t &irow_recv);

    void annihilate();
};


#endif //M7_SOLVER_H
