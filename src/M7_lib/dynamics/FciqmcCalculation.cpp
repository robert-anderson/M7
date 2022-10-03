//
// Created by Robert John Anderson on 2020-02-11.
//

#include <M7_lib/io/Logging.h>
#include <M7_lib/field/Mbf.h>

#include "FciqmcCalculation.h"
#include "M7_lib/propagator/Propagators.h"

FciqmcCalculation::FciqmcCalculation(const conf::Document &opts) :
        m_opts(opts), m_ham({opts.m_hamiltonian, opts.m_basis}),
        m_wf(opts, {m_ham.m_basis, m_ham.default_particles(opts.m_particles)}),
        m_prop(props::get(m_ham, opts, m_wf)) {
    buffered::Mbf ref_mbf(m_wf.m_sector);

    mbf::set(ref_mbf, m_wf.m_sector.particles(), opts.m_reference.m_mbf_init, 0ul);

    auto ref_energy = m_ham.get_energy(ref_mbf);
    TableBase::Loc ref_loc = {m_wf.m_dist.irank(ref_mbf), 0ul};
    m_wf.create_row(0, ref_mbf, ref_energy, v_t<bool>(m_wf.npart(), true));
    if (ref_loc.is_mine()){
        auto ref_walker = m_wf.m_store.m_row;
        ref_walker.jump(ref_loc.m_irec);
        for (uint_t ipart=0ul; ipart<m_wf.npart(); ++ipart)
            m_wf.set_weight(ref_walker, ipart, opts.m_wavefunction.m_nw_init);
    }
    m_prop->m_shift.m_values = ref_energy+opts.m_shift.m_init;
    Solver solver(opts, *m_prop, m_wf, ref_loc);
    solver.execute(opts.m_propagator.m_ncycle);
}
