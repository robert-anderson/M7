//
// Created by Robert John Anderson on 2020-02-11.
//

#ifndef M7_FCIQMCCALCULATION_H
#define M7_FCIQMCCALCULATION_H


#include <omp.h>
#include <src/core/util/Timer.h>
#include "src/core/hamiltonian/AbInitioHamiltonian.h"
#include "src/core/io/FciqmcStatsFile.h"
#include "src/core/io/Options.h"
#include "src/core/fermion/DecodedDeterminant.h"
#include "Propagator.h"
#include "Wavefunction.h"

class FciqmcCalculation {
public:
    const Options m_input;
    RankAllocator<DeterminantElement> m_rank_allocator;
    std::unique_ptr<FciqmcStatsFile> m_stats_file = nullptr;
    std::unique_ptr<Hamiltonian> m_ham;
    Determinant m_reference;
    std::unique_ptr<Propagator> m_prop;
    Wavefunction m_wf;
    Timer m_timer;

    explicit FciqmcCalculation(const Options &input);

    void execute();

    void write_iter_stats(size_t icycle);

};


#endif //M7_FCIQMCCALCULATION_H
