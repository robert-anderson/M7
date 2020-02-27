//
// Created by Robert John Anderson on 2020-02-11.
//

#include "FciqmcCalculation.h"
#include "src/io/Logging.h"


FciqmcCalculation::FciqmcCalculation(const InputOptions &input) :
        m_input(input), m_rank_allocator(input.nload_balance_block), m_stats_file() {
    m_ham = std::make_unique<AbInitioHamiltonian>(defs::assets_root + "/DHF_Be_STO-3G/FCIDUMP");
    m_prop = std::make_unique<ExactPropagator>(m_ham, m_rank_allocator, input.tau_initial, input.shift_initial);
    m_psi = std::make_unique<Wavefunction>(
            m_ham->guess_reference(input.spin_level),
            m_input.nwalker_initial,
            (size_t) (m_input.nwalker_target * m_input.walker_factor_initial),
            (size_t) (m_input.nwalker_target * m_input.buffer_factor_initial),
            (size_t) (m_input.nwalker_target * m_input.buffer_factor_initial * mpi::nrank()));
}

void FciqmcCalculation::execute(size_t ncycle){
    logger::write("Starting FCIQMC main loop.");
    for(size_t icycle=0ul; icycle<ncycle; ++icycle){
        m_psi->propagate(m_prop);
        m_psi->communicate();
        m_psi->consolidate_incoming_weight();
        m_psi->annihilate();
        m_prop->update(icycle, m_psi->m_norm_growth_rate);
        m_psi->write_iter_stats();
    }
}