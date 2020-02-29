//
// Created by rja on 28/02/2020.
//

#include "FciqmcStatsFile.h"

FciqmcStatsFile::FciqmcStatsFile(const InputOptions &input) : StatsFile(input.stats_path) {
    /*
     * Setup stats columns
     */
    m_cycle_number =
            add_column<size_t>("Cycle number");

    m_diagonal_shift =
            add_column<defs::ham_comp_t>("Diagonal shift");

    m_timestep =
            add_column<double>("Timestep");

    m_reference_projected_energy_numerator =
            add_column<defs::ham_t>("Reference projected energy numerator");

    m_reference_weight =
            add_column<defs::ham_t>("Reference weight");

    m_reference_energy =
            add_column<defs::ham_t>("Reference energy");

    m_wavefunction_l2_norm =
            add_column<defs::ham_comp_t>("Wavefunction L2 norm");

    write_header();
}
