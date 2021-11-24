//
// Created by Robert John Anderson on 2020-02-11.
//

#include "FciqmcCalculation.h"
#include "src/core/io/Logging.h"
#include "ExactPropagator.h"
#include "StochasticPropagator.h"
#include "Propagators.h"

FciqmcCalculation::FciqmcCalculation(const fciqmc_config::Document &opts) :
        m_opts(opts), m_ham(opts.m_hamiltonian), m_wf(opts, m_ham.m_bd),
        m_prop(props::get(m_ham, opts, m_wf.m_format)) {
    buffered::Mbf ref_onv(m_ham.m_bd);
    //m_ham.set_mbf(ref_onv, opts.m_reference.m_init_mbf_neel, opts.m_wavefunction.m_spin_restrict);
    ref_onv = {0, 0, 0, 5, 0};
    auto ref_energy = m_ham.get_energy(ref_onv);
    TableBase::Loc ref_loc = {m_wf.get_rank(ref_onv), 0ul};
    m_wf.create_row(0, ref_onv, ref_energy, std::vector<bool>(m_wf.npart(), true));
    if (ref_loc.is_mine()){
        for (size_t ipart=0ul; ipart<m_wf.npart(); ++ipart)
            m_wf.set_weight(ipart, opts.m_wavefunction.m_nw_init.get());
    }
    m_prop->m_shift.m_values = ref_energy+opts.m_shift.m_init;
    Solver solver(opts, *m_prop, m_wf, ref_loc);
    solver.execute(opts.m_propagator.m_ncycle);
}