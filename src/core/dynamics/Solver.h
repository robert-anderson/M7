//
// Created by rja on 10/11/2020.
//

#ifndef M7_SOLVER_H
#define M7_SOLVER_H

#include <src/core/hamiltonian/Hamiltonian.h>
#include <src/core/util/Timer.h>
#include "Reference.h"
#include "src/core/table/Communicator.h"
#include "src/core/field/Fields.h"
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

    Solver(Propagator &prop, Wavefunction &wf, Table::Loc ref_loc);

    void execute(size_t niter=1);

    void begin_cycle();

    void propagate_row(const size_t& irow);

    void loop_over_occupied_onvs();

    void annihilate_row(const size_t &irow_recv);

    void loop_over_spawned();

    void end_cycle();

    void output_stats();
};




//    mpi::barrier();
//    //m_propagation_timer.pause();
//
//    /*
//     * the numerator of the reference projected energy Rayleigh quotient is treated like a
//     * delta variable, but is an instantaneous principal variable determined by the pre-
//     * propagation walker distribution. This also enables the instantaneous Rayleigh quotient
//     * to be computed. This is a pure convenience to aid a "quick look" at the energy estimator
//     * in the stats file - the numerator and denominator should be independently analyzed to
//     * obtain the estimate. All of this is handled in the Reference class, so defer to the
//     * method implemented therein.
//     */
//    m_reference.synchronize();
//
//#ifdef VERBOSE_DEBUGGING
//    std::cout << consts::verb << consts::chevs << "END OF PROPAGATION LOOP CHECKS" << std::endl;
//std::cout << consts::verb << "free rows found in walker list:    " << nrow_free << std::endl;
//std::cout << consts::verb << "free rows after propagation:       " << m_data.nzero_rows() << std::endl;
//std::cout << consts::verb << "occupied determinants before loop: " << m_nocc_det.local() << std::endl;
//std::cout << consts::verb << "delta in occupied determinants:    " << m_nocc_det.m_delta.local() << std::endl;
//std::cout << consts::verb << "map size:                          " << m_data.map_size() << std::endl;
//std::cout << consts::verb << "high water mark:                   " << m_data.high_water_mark(0) << std::endl;
//#endif
//
//#if 0
//    #ifndef NDEBUG
//        ASSERT(nrow_free - m_nocc_det.m_delta.local() == m_data.nzero_rows())
//
//        auto chk_hwm = nrow_free + m_nocc_det.local();
//        if (chk_hwm != m_data.high_water_mark(0)) {
//            m_data.print();
//            auto chk_nrow_in_free_stack = m_data.nrow_in_free_stack();
//            std::cout << "free rows in walker list " << chk_nrow_in_free_stack << std::endl;
//        }
//        ASSERT(chk_hwm == m_data.high_water_mark(0))
//        ASSERT(chk_hwm - m_data.nzero_rows() == m_data.map_size())
//        ASSERT(consts::floats_equal(chk_proj_energy, m_reference.proj_energy()))
//#endif
//#endif
//
//
//}

#endif //M7_SOLVER_H
