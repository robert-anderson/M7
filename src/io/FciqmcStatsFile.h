//
// Created by rja on 28/02/2020.
//

#ifndef M7_FCIQMCSTATSFILE_H
#define M7_FCIQMCSTATSFILE_H

#include "StatsFile.h"
#include "InputOptions.h"

struct FciqmcStatsFile : public StatsFile {
public:
    StatColumn *m_cycle_number;
    StatColumn *m_diagonal_shift;
    StatColumn *m_timestep;
    StatColumn *m_reference_projected_energy_numerator;
    StatColumn *m_reference_weight;
    StatColumn *m_reference_energy;
    StatColumn *m_wavefunction_l2_norm;
    StatColumn *m_ninitiator;
    StatColumn *m_aborted_weight;

    explicit FciqmcStatsFile(const InputOptions &input);
};


#endif //M7_FCIQMCSTATSFILE_H
