//
// Created by Robert John Anderson on 2020-02-11.
//

#include "FciqmcCalculation.h"
#include "src/core/io/Logging.h"

FciqmcCalculation::FciqmcCalculation(const InputOptions &input) :
        m_input(input), m_rank_allocator(input.nload_balance_block),
        m_stats_file(input)
        {
    m_stats_file.m_ninitiator() = 13;
    m_stats_file.m_ref_weight() = {1,2};
    m_stats_file.flush();
    //m_ham = std::make_unique<AbInitioHamiltonian>(input.fcidump_path);
    //m_prop = std::make_unique<ExactPropagator>(input, m_ham, m_rank_allocator);
    //m_prop = std::make_unique<StochasticPropagator>(input, m_ham, m_rank_allocator);
    //auto reference = m_ham->guess_reference(input.spin_restrict);
    //m_prop->m_shift += m_ham->get_energy(reference);
    //m_psi = std::make_unique<Wavefunction>(input, m_prop, reference);
}

void FciqmcCalculation::execute(){
    logger::write("Starting FCIQMC main loop.");
    for(size_t icycle=0ul; icycle<m_input.ncycle; ++icycle){
        /*
        logger::write("iteration "+std::to_string(icycle));
        logger::write("\npropagating...");
        m_psi->propagate(m_prop);
        logger::write("\ncommunicating...");
        m_psi->communicate();
        m_psi->consolidate_incoming_weight();
        logger::write("\nannihilating...");
        m_psi->annihilate(m_prop);
        m_prop->update(icycle, m_psi->norm(), m_psi->m_norm_growth_rate);
        write_iter_stats(icycle);
         */
    }
}

void FciqmcCalculation::write_iter_stats(size_t icycle) {
    /*
    m_stats_file.m_cycle_number->write(icycle);
    m_prop->write_iter_stats(m_stats_file);
    m_psi->write_iter_stats(m_stats_file);
    m_stats_file.flush();
     */
}
