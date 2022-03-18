//
// Created by Robert John Anderson on 2020-02-11.
//

#ifndef M7_FCIQMCCALCULATION_H
#define M7_FCIQMCCALCULATION_H


#include <config/FciqmcConfig.h>
#include <wavefunction/Wavefunction.h>

#include "StochasticPropagator.h"
#include "Solver.h"

class FciqmcCalculation {
public:
    const fciqmc_config::Document& m_opts;
    Hamiltonian m_ham;
    Wavefunction m_wf;
    std::unique_ptr<Propagator> m_prop;

    explicit FciqmcCalculation(const fciqmc_config::Document& opts);
};



#if 0

#include <omp.h>
#include <util/Timer.h>
#include "hamiltonian/AbInitioHamiltonian.h"
#include "io/FciqmcStatsFile.h"
#include "io/ParallelizationStatsFile.h"
#include "io/Options.h"
#include "basis/DecodedDeterminant.h"
#include "Propagator.h"
#include "Wavefunction.h"

class FciqmcCalculation {
public:
    Epoch m_vary_shift;
    Epoch m_semi_stochastic;
    const Options m_input;
    RankAllocator<DeterminantElement> m_rank_allocator;
    std::unique_ptr<FciqmcStatsFile> m_stats_file = nullptr;
    std::unique_ptr<ParallelizationStatsFile> m_parallel_stats_file = nullptr;

    std::unique_ptr<FermionHamiltonian> m_ham;
    FermionOnv m_reference;
    std::unique_ptr<Propagator> m_prop;
    Wavefunction m_wf;
    Timer m_timer;

    explicit FciqmcCalculation(const Options &input);

    void execute();

    void write_iter_stats(size_t icycle);

    FermionOnv initial_reference(const std::unique_ptr<FermionHamiltonian>& ham, const Options& input){
        if(m_input.initial_reference_det.empty()) {
            return m_ham->guess_reference(input.spin_restrict);
        }
        else {
            FermionOnv result(ham->nsite());
            result.set(input.initial_reference_det);
            return result;
        }
    }

};


#endif //M7_FCIQMCCALCULATION_H
#endif //M7_FCIQMCCALCULATION_H
