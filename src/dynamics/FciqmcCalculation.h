//
// Created by Robert John Anderson on 2020-02-11.
//

#ifndef M7_FCIQMCCALCULATION_H
#define M7_FCIQMCCALCULATION_H


#include <omp.h>
#include <src/io/StatsFile.h>
#include <src/hamiltonian/AbInitioHamiltonian.h>
#include "src/io/InputOptions.h"
#include "src/hamiltonian/Hamiltonian.h"
#include "Propagator.h"
#include "Wavefunction.h"
#include "ExactPropagator.h"

class FciqmcCalculation {
    const InputOptions m_input;
    RankAllocator<Determinant> m_rank_allocator;
    std::unique_ptr<Hamiltonian> m_ham;
    std::unique_ptr<Propagator> m_prop;
    std::unique_ptr<Wavefunction> m_psi;
    StatsFile m_stats_file;
public:
    explicit FciqmcCalculation(const InputOptions &input);

    void execute(size_t ncycle);

};


#endif //M7_FCIQMCCALCULATION_H
