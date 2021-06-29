//
// Created by Robert John Anderson on 2020-02-11.
//

#include "FciqmcCalculation.h"
#include "src/core/io/Logging.h"
#include "ExactPropagator.h"
#include "StochasticPropagator.h"
#include "Propagators.h"

FciqmcCalculation::FciqmcCalculation(const fciqmc_config::Document& opts) :
        m_opts(opts), m_ham(opts.m_hamiltonian), m_wf(opts, m_ham.nsite()),
        m_prop(props::get(m_ham, opts, m_wf.m_format))  {
    buffered::Onv<> ref_onv(m_ham.nsite());
    m_ham.set_hf_onv(ref_onv, opts.m_wavefunction.m_spin_restrict);
    auto ref_energy = m_ham.get_energy(ref_onv);
    TableBase::Loc ref_loc = {m_wf.get_rank(ref_onv), 0ul};
    if (ref_loc.is_mine()) {
        m_wf.create_row(0, ref_onv, ref_energy, std::vector<bool>(m_wf.npart(), true));
        m_wf.set_weight(0, ref_energy);
    }
    m_prop->m_shift.m_values = ref_energy;
    Solver solver(opts, *m_prop, m_wf, ref_loc);
    solver.execute(opts.m_propagator.m_ncycle);
}

#if 0

FciqmcCalculation::FciqmcCalculation(const Options &input) :
        m_input(input), m_vary_shift("variable shift"), m_semi_stochastic("semi-stochastic"),
        m_rank_allocator(input.nload_balance_block_per_rank*mpi::nrank(),
            input.load_balance_period, &m_vary_shift),
        m_ham(std::unique_ptr<AbInitioHamiltonian>(
                new AbInitioHamiltonian(input.fcidump_path, input.fcidump_spin_major))),
        m_reference(initial_reference(m_ham, input)),
        m_wf(this) {

    logger::write("Initializing FCIQMC Calculation...");
    if (mpi::i_am_root()) m_stats_file = std::unique_ptr<FciqmcStatsFile>(new FciqmcStatsFile(m_input));
    m_parallel_stats_file = std::unique_ptr<ParallelizationStatsFile>(new ParallelizationStatsFile(m_input));

    if (input.exact_propagation) {
        m_prop = std::unique_ptr<ExactPropagator>(new ExactPropagator(this));
    } else {
        m_prop = std::unique_ptr<StochasticPropagator>(new StochasticPropagator(this));
    }
    m_prop->m_values += m_ham->get_energy(m_reference);

    logger::write("Distributed memory parallelization: " + std::to_string(mpi::nrank()) + " MPI ranks");
    logger::write("Reference determinant was detected to be: " + m_reference.to_string());
}

void FciqmcCalculation::execute() {
    logger::write("Starting FCIQMC main loop.");
    for (size_t icycle = 0ul; icycle < m_input.ncycle; ++icycle) {
        mpi::barrier();
        m_timer.unpause();
        //std::cout << "iteration " + std::to_string(icycle) << std::endl;
        m_wf.update(icycle);
        //std::cout << "\npropagating..." << std::endl;
        m_wf.propagate();
        //std::cout << "\ncommunicating..." << std::endl;
        m_wf.communicate();
        m_wf.consolidate_incoming_weight();
        //std::cout << "\nannihilating..." << std::endl;
        m_wf.annihilate();
        m_wf.synchronize();
        //std::cout << "\nupdating propagator..." << std::endl;
        m_prop->update(icycle, m_wf.m_nwalker.reduced(), m_wf.m_nwalker_growth_rate);
        mpi::barrier();
        m_timer.pause();
        //std::cout << "\nwriting stats..." << std::endl;
        write_iter_stats(icycle);
        m_wf.write_iter_stats(m_stats_file.get());
    }
}

void FciqmcCalculation::write_iter_stats(size_t icycle) {
    m_parallel_stats_file->m_cycle_number.write(icycle);
    m_parallel_stats_file->m_synchronization_overhead.write(0.0);
    m_parallel_stats_file->flush();
    if (!mpi::i_am_root()) return;
    m_stats_file->m_cycle_number.write(icycle);
    m_stats_file->m_iter_time.write(m_timer.lap());
    m_prop->write_iter_stats(m_stats_file.get());
    m_wf.write_iter_stats(m_stats_file.get());
    m_stats_file->flush();
}

#endif
