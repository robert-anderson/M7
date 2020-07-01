//
// Created by Robert John Anderson on 2020-02-11.
//

#include "FciqmcCalculation.h"
#include "src/core/io/Logging.h"
#include "ExactPropagator.h"
#include "StochasticPropagator.h"


FciqmcCalculation::FciqmcCalculation(const Options &input) :
    m_input(input), m_rank_allocator(input.nload_balance_block),
    m_ham(std::unique_ptr<AbInitioHamiltonian>(new AbInitioHamiltonian(input.fcidump_path))),
    m_reference(m_ham->guess_reference(input.spin_restrict)),
    m_wf(this), m_scratch(std::unique_ptr<FciqmcScratch>(new FciqmcScratch(m_reference))) {

    if(mpi::i_am_root()) m_stats_file = std::unique_ptr<FciqmcStatsFile>(new FciqmcStatsFile(m_input));

    if (input.exact_propagation) {
        m_prop = std::unique_ptr<ExactPropagator>(new ExactPropagator(this));
    } else {
        m_prop = std::unique_ptr<StochasticPropagator>(new StochasticPropagator(this));
    }
    m_prop->m_shift += m_ham->get_energy(m_reference);

    logger::write("Initializing FCIQMC Calculation...");
    logger::write("Shared memory parallelization: "+std::to_string(omp_get_max_threads())+" OpenMP threads");
    logger::write("Distributed memory parallelization: "+std::to_string(mpi::nrank())+" MPI ranks");
    logger::write("Reference determinant was detected to be: " + m_reference.to_string());
}

void FciqmcCalculation::execute() {
    logger::write("Starting FCIQMC main loop.");
    for (size_t icycle = 0ul; icycle < m_input.ncycle; ++icycle) {
        //std::cout << "iteration " + std::to_string(icycle) << std::endl;
        m_wf.update(icycle);
        //std::cout << "\npropagating..." << std::endl;
        m_wf.propagate();
        //std::cout << "\ncommunicating..." << std::endl;
        m_wf.communicate();
        m_wf.consolidate_incoming_weight();
        //std::cout << "\nannihilating..." << std::endl;
        m_wf.annihilate();
        //std::cout << "\nupdating propagator..." << std::endl;
        m_prop->update(icycle, m_wf.m_nw.reduced(), m_wf.m_nw_growth_rate);
        //std::cout << "\nwriting stats..." << std::endl;
        write_iter_stats(icycle);
    }
}

void FciqmcCalculation::write_iter_stats(size_t icycle) {
    if (!mpi::i_am_root()) return;
    m_stats_file->m_cycle_number.write(icycle);
    m_prop->write_iter_stats(m_stats_file.get());
    m_wf.write_iter_stats(m_stats_file.get());
    m_stats_file->flush();
}
