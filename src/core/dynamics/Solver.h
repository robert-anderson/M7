//
// Created by rja on 10/11/2020.
//

#ifndef M7_SOLVER_H
#define M7_SOLVER_H

#include <src/core/hamiltonian/Hamiltonian.h>
#include <src/core/util/Timer.h>
#include "Reference.h"
#include "src/core/table/Communicator.h"
#include "src/core/io/FciqmcStatsFile.h"
#include "src/core/io/ParallelStatsFile.h"
#include "Propagator.h"

class Solver {

    size_t m_icycle = 0ul;
    Propagator &m_prop;
    const Options &m_opts;
    Wavefunction &m_wf;
    Reference m_reference;

    StatsFile<FciqmcStatsSpecifier>::ptr_t m_stats;
    StatsFile<ParallelStatsSpecifier>::ptr_t m_parallel_stats;

    /*
     * Timers for the main parts of the solver
     */
    // whole cycle
    Timer m_cycle_timer;
    // whole loop over occupied rows
    Timer m_propagate_timer;
    // individual iterations over an occupied row
    Timer m_spawning_timer;
    // time waited at MPI_Barrier
    Timer m_synchronization_timer;
    // time taken to complete communication of spawned walkers
    Timer m_communicate_timer;
    // time taken to complete whole annihilation loop
    Timer m_annihilate_timer;

    /*
     * Sanity checking variables
     */
    defs::wf_t m_chk_nwalker_local = 0.0;
    size_t m_chk_ninitiator_local = 0ul;

public:

    Solver(Propagator &prop, Wavefunction &wf, TableBase::Loc ref_loc);

    void execute(size_t niter=1);

    void begin_cycle();

    void propagate_row();

    void loop_over_occupied_onvs();

    void annihilate_row();

    void loop_over_spawned();

    void end_cycle();

    void output_stats();
};

#endif //M7_SOLVER_H
