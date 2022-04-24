//
// Created by Robert John Anderson on 2020-02-11.
//

#include <M7_lib/io/Logging.h>
#include <M7_lib/field/Mbf.h>

#include "FciqmcCalculation.h"
#include "Propagators.h"

FciqmcCalculation::FciqmcCalculation(const fciqmc_config::Document &opts) :
        m_opts(opts), m_ham(opts.m_hamiltonian), m_wf(opts, m_ham.m_hs),
        m_prop(props::get(m_ham, opts, m_wf.m_format)) {
    buffered::Mbf ref_mbf(m_ham.m_hs);

    mbf::set(ref_mbf, opts.m_reference.m_mbf_init, 0ul);

    auto ref_energy = m_ham.get_energy(ref_mbf);
    TableBase::Loc ref_loc = {m_wf.get_rank(ref_mbf), 0ul};
    m_wf.create_row(0, ref_mbf, ref_energy, std::vector<bool>(m_wf.npart(), true));
    if (ref_loc.is_mine()){
        for (size_t ipart=0ul; ipart<m_wf.npart(); ++ipart)
            m_wf.set_weight(ipart, opts.m_wavefunction.m_nw_init.get());
    }
    m_prop->m_shift.m_values = ref_energy+opts.m_shift.m_init;
    Solver solver(opts, *m_prop, m_wf, ref_loc);
    solver.execute(opts.m_propagator.m_ncycle);
}
