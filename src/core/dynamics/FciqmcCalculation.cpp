//
// Created by Robert John Anderson on 2020-02-11.
//

#include "FciqmcCalculation.h"
#include "src/core/io/Logging.h"
#include "ExactPropagator.h"
#include "StochasticPropagator.h"


FciqmcCalculation::FciqmcCalculation(const Options &opts) :
        m_opts(opts), m_ham(opts), m_prop(m_ham, opts), m_wf(opts, m_ham.nsite(), &m_prop.m_variable_shift) {
    m_wf.expand(size_t(opts.walker_factor_initial * opts.nwalker_target),
                size_t(opts.buffer_factor_initial * opts.nwalker_target));
    auto ref_det = m_ham.guess_reference(opts.spin_restrict);
    auto ref_energy = m_ham.get_energy(ref_det);
    m_prop.m_shift = ref_energy;//benchmark;
    Solver solver(m_prop, m_wf, ref_det);
    for (size_t i = 0ul; i < opts.ncycle; ++i) {
        solver.execute();
    }
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
    m_prop->m_shift += m_ham->get_energy(m_reference);

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
    m_parallel_stats_file->m_synchronization_wait_time.write(0.0);
    m_parallel_stats_file->flush();
    if (!mpi::i_am_root()) return;
    m_stats_file->m_cycle_number.write(icycle);
    m_stats_file->m_iter_time.write(m_timer.lap());
    m_prop->write_iter_stats(m_stats_file.get());
    m_wf.write_iter_stats(m_stats_file.get());
    m_stats_file->flush();
}

#endif
