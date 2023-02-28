//
// Created by Robert John Anderson on 2020-02-11.
//

#include <M7_lib/io/Logging.h>
#include <M7_lib/field/Mbf.h>

#include "FciqmcCalculation.h"
#include "M7_lib/propagator/Propagators.h"

void fciqmc::run(const conf::Document &opts) {
    /*
     * sum of weighted many-body operator products determining the energies and transition amplitudes between MBFs
     */
    const Hamiltonian ham({opts.m_hamiltonian, opts.m_basis, opts.m_particles});
    /*
     * distributed solution vectors
     */
    wf::Fci wf(opts, {ham.m_basis, ham.default_particles(opts.m_particles)});
    /*
     * propagates the system, either exactly or stochastically
     */
    std::unique_ptr<Propagator> m_prop = props::get(ham, opts, wf);
    /*
     * create the reference MBF and add it to the walker table
     */
    buffered::Mbf ref_mbf(wf.m_sector);
    mbf::set(ref_mbf, wf.m_sector.particles(), opts.m_reference.m_mbf_init, 0ul);

    auto ref_energy = ham.get_energy(ref_mbf);
    TableBase::Loc ref_loc = {wf.m_dist.irank(ref_mbf), 0ul};
    wf.create_row(0, ref_mbf, ref_energy, v_t<bool>(wf.npart(), true));
    if (ref_loc.is_mine()) {
        auto ref_walker = wf.m_store.m_row;
        ref_walker.jump(ref_loc.m_irec);
        for (uint_t ipart = 0ul; ipart < wf.npart(); ++ipart)
            wf.set_weight(ref_walker, ipart, wf_t(opts.m_wavefunction.m_nw_init));
    }
    m_prop->m_shift.m_values = ref_energy + opts.m_shift.m_init;
    wf.m_nwalker.all_sum();
    Solver solver(opts, *m_prop, wf, ref_loc);
    solver.execute(opts.m_propagator.m_ncycle);
}
