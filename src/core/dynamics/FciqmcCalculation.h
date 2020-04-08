//
// Created by Robert John Anderson on 2020-02-11.
//

#ifndef M7_FCIQMCCALCULATION_H
#define M7_FCIQMCCALCULATION_H


#include <omp.h>
#include "src/core/hamiltonian/AbInitioHamiltonian.h"
#include "src/core/io/FciqmcStatsFile.h"
#include "src/core/io/InputOptions.h"
#include "Propagator.h"
#include "Wavefunction.h"

class FciqmcCalculation {
    const InputOptions m_input;
    RankAllocator<DeterminantElement> m_rank_allocator;
    FciqmcStatsFile m_stats_file;

#if 0

    std::unique_ptr<Hamiltonian> m_ham;
    std::unique_ptr<Propagator> m_prop;
    std::unique_ptr<Wavefunction> m_psi;
#endif //M7_FCIQMCCALCULATION_H
public:
    explicit FciqmcCalculation(const InputOptions &input);

    void execute();

    void write_iter_stats(size_t icycle);

};


#endif //M7_FCIQMCCALCULATION_H
