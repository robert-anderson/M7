//
// Created by Robert John Anderson on 2020-02-04.
//

#ifndef M7_WAVEFUNCTION_H
#define M7_WAVEFUNCTION_H

#include <src/core/io/StatsFile.h>
#include <src/core/io/FciqmcStatsFile.h>
#include "WalkerList.h"
#include "src/heatbath/HeatBathSampler.h"
#include "WalkerCommunicator.h"
#include "Propagator.h"
#include "RankAllocator.h"

#if 0
class Wavefunction {
    WalkerList m_walker_list;
    //WalkerCommunicator m_walker_communicator;

    Determinant m_reference;
    defs::ham_t m_reference_energy_numerator;
    defs::ham_t m_reference_weight;
    defs::ham_comp_t m_aborted_weight;
    size_t m_ninitiator = 0;

public:

    const InputOptions &m_input;
    defs::ham_comp_t m_square_norm;
    defs::ham_comp_t m_delta_square_norm;
    size_t m_noccupied_determinant;
    double m_norm_growth_rate = 0;

    Wavefunction(const InputOptions &input, const std::unique_ptr<Propagator> &propagator,
                 const Determinant &reference);

    void propagate(std::unique_ptr<Propagator> &propagator);

    void communicate();

    void consolidate_incoming_weight() {
        // TODO
    }

    void annihilate(const std::unique_ptr<Propagator> &propagator);

    void write_iter_stats(FciqmcStatsFile &stats_file);

    defs::ham_comp_t norm() const;

};


#endif //M7_WAVEFUNCTION_H
#endif //M7_WAVEFUNCTION_H
