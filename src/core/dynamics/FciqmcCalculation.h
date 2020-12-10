//
// Created by Robert John Anderson on 2020-02-11.
//

#ifndef M7_FCIQMCCALCULATION_H
#define M7_FCIQMCCALCULATION_H


#include <src/core/io/Options.h>
#include "Wavefunction.h"
#include "StochasticPropagator.h"

//class FciqmcCalculation {
//public:
//    const Options m_opts;
//    Hamiltonian<defs::enable_bosons> m_ham;
//    Wavefunction m_wf;
//    StochasticPropagator m_prop;
//
//    FciqmcCalculation(const Options& opts):m_opts(opts),
//    m_ham(opts), m_wf(opts, m_ham.nsite()), m_prop(m_ham, opts){
//        m_wf.expand(10, 800);
//        auto ref_energy = m_ham.get_energy(fonv);
//        prop.m_shift = ref_energy;//benchmark;
//
//    }
//    Solver solver(prop, wf, fonv);
//
//    std::cout << "Reference Energy: " << ref_energy << std::endl;
//
//    for (size_t i = 0ul; i < 10000; ++i) {
//        solver.execute();
//        std::cout << i << " " << wf.m_walkers.m_hwm << " " << std::sqrt(wf.square_norm()) << std::endl;
//    }
//};



#if 0

#include <omp.h>
#include <src/core/util/Timer.h>
#include "src/core/hamiltonian/AbInitioHamiltonian.h"
#include "src/core/io/FciqmcStatsFile.h"
#include "src/core/io/ParallelizationStatsFile.h"
#include "src/core/io/Options.h"
#include "src/core/basis/DecodedDeterminant.h"
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
