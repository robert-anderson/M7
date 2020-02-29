//
// Created by Robert John Anderson on 2020-02-11.
//

#include "FciqmcCalculation.h"
#include "src/io/Logging.h"


FciqmcCalculation::FciqmcCalculation(const InputOptions &input) :
        m_input(input), m_rank_allocator(input.nload_balance_block),
        m_stats_file(input)
        {
    m_ham = std::make_unique<AbInitioHamiltonian>(defs::assets_root + "/DHF_Be_STO-3G/FCIDUMP");
    m_prop = std::make_unique<ExactPropagator>(input, m_ham, m_rank_allocator);
    auto reference = m_ham->guess_reference(input.spin_level);
    m_prop->m_shift += m_ham->get_energy(reference);
    m_psi = std::make_unique<Wavefunction>(
            m_prop,
            reference,
            m_input.wf_norm_initial,
            (size_t) (m_input.wf_norm_target * m_input.walker_factor_initial),
            (size_t) (m_input.wf_norm_target * m_input.buffer_factor_initial),
            (size_t) (m_input.wf_norm_target * m_input.buffer_factor_initial * mpi::nrank()));
}

void FciqmcCalculation::execute(size_t ncycle){
    logger::write("Starting FCIQMC main loop.");
    for(size_t icycle=0ul; icycle<ncycle; ++icycle){
        m_psi->propagate(m_prop);
        m_psi->communicate();
        m_psi->consolidate_incoming_weight();
        m_psi->annihilate(m_prop);
        m_prop->update(icycle, m_psi->norm(), m_psi->m_norm_growth_rate);
        write_iter_stats(icycle);
    }
}

void FciqmcCalculation::write_iter_stats(size_t icycle) {
    m_stats_file.m_cycle_number->write(icycle);
    m_prop->write_iter_stats(m_stats_file);
    m_psi->write_iter_stats(m_stats_file);
    m_stats_file.flush();
}
