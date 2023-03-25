//
// Created by Robert John Anderson on 2020-02-11.
//

#include <M7_lib/io/Logging.h>
#include <M7_lib/field/Mbf.h>

#include "FciqmcCalculation.h"
#include "M7_lib/propagator/Propagators.h"

void fciqmc::run(const conf::Document &opts) {
    mpi::barrier();
    /*
     * sum of weighted many-body operator products determining the energies and transition amplitudes between MBFs
     */
    const Hamiltonian ham({opts.m_hamiltonian, opts.m_basis, opts.m_particles});
    /*
     * distributed solution vectors
     */
    wf::Vectors wf(opts, ham);
    /*
     * propagates the system, either exactly or stochastically
     */
    std::unique_ptr<Propagator> m_prop = props::get(ham, opts, wf);

    const ham_comp_t ground_ref_energy = wf.m_refs[0].hdiag();
    m_prop->m_shift.m_values = ground_ref_energy + opts.m_shift.m_init;
    wf.m_stats.m_nwalker.all_sum();
    Solver solver(opts, *m_prop, wf);
    solver.execute(opts.m_propagator.m_ncycle);
}
