//
// Created by Robert John Anderson on 2020-02-11.
//

#include "FciqmcCalculation.h"
#include "src/core/io/Logging.h"
#include "ExactPropagator.h"


FciqmcCalculation::FciqmcCalculation(const InputOptions &input) :
    m_input(input), m_rank_allocator(input.nload_balance_block),
    m_stats_file(input),
    m_ham(std::make_unique<AbInitioHamiltonian>(input.fcidump_path)),
    m_reference(m_ham->guess_reference(input.spin_restrict)),
    m_wf(this), m_scratch(std::make_unique<FciqmcScratch>(m_reference)) {
    m_prop = std::make_unique<ExactPropagator>(this),
        //m_prop = std::make_unique<StochasticPropagator>(input, m_ham, m_rank_allocator);
        m_prop->m_shift += m_ham->get_energy(m_reference);

    logger::write("Initializing FCIQMC Calculation...");
    logger::write("Reference determinant was detected to be: " + m_reference.to_string());

}

void FciqmcCalculation::execute() {
    logger::write("Starting FCIQMC main loop.");
    for (size_t icycle = 0ul; icycle < m_input.ncycle; ++icycle) {
        logger::write("iteration " + std::to_string(icycle));
        //logger::write("\npropagating...");
        m_wf.propagate();
        //logger::write("\ncommunicating...");
        m_wf.communicate();
        m_wf.consolidate_incoming_weight();
        //logger::write("\nannihilating...");
        m_wf.annihilate();
        m_prop->update(icycle, m_wf.m_nw, m_wf.m_nw_growth_rate);
        write_iter_stats(icycle);
    }
}

void FciqmcCalculation::write_iter_stats(size_t icycle) {
    m_stats_file.m_cycle_number() = icycle;
    m_prop->write_iter_stats(m_stats_file);
    m_wf.write_iter_stats(m_stats_file);
    m_stats_file.flush();
}
