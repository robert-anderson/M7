//
// Created by Robert John Anderson on 2020-02-11.
//

#ifndef M7_FCIQMCCALCULATION_H
#define M7_FCIQMCCALCULATION_H


#include <omp.h>
#include <src/core/util/Timer.h>
#include "src/core/hamiltonian/AbInitioHamiltonian.h"
#include "src/core/io/FciqmcStatsFile.h"
#include "src/core/io/ParallelizationStatsFile.h"
#include "src/core/io/Options.h"
#include "src/core/fermion/DecodedDeterminant.h"
#include "Propagator.h"
#include "Wavefunction.h"

class FciqmcCalculation {
public:
    const Options m_input;
    Epoch m_vary_shift;
    Epoch m_semi_stochastic;
    RankAllocator<DeterminantElement> m_rank_allocator;
    std::unique_ptr<FciqmcStatsFile> m_stats_file = nullptr;
    std::unique_ptr<ParallelizationStatsFile> m_parallel_stats_file = nullptr;

    std::unique_ptr<Hamiltonian> m_ham;
    Determinant m_reference;
    std::unique_ptr<Propagator> m_prop;
    Wavefunction m_wf;
    Timer m_timer;

    explicit FciqmcCalculation(const Options &input);

    void execute();

    void write_iter_stats(size_t icycle);

    Determinant initial_reference(const std::unique_ptr<Hamiltonian>& ham, const Options& input){
        if(m_input.initial_reference_det.empty()) {
            return m_ham->guess_reference(input.spin_restrict);
        }
        else {
            Determinant result(ham->nsite());
            result.set(input.initial_reference_det);
            return result;
        }
    }

};


#endif //M7_FCIQMCCALCULATION_H
