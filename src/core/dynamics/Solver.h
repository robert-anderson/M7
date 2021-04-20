//
// Created by rja on 10/11/2020.
//

#ifndef M7_SOLVER_H
#define M7_SOLVER_H

#include <src/core/hamiltonian/Hamiltonian.h>
#include <src/core/util/Timer.h>
#include <src/core/observables/AverageCoefficients.h>
#include <src/core/observables/UniformTwf.h>
#include <src/core/observables/StaticTwf.h>
#include <src/core/observables/MevGroup.h>
#include <src/core/io/FciqmcStats.h>
#include "Reference.h"
#include "src/core/table/Communicator.h"
#include "src/core/io/FciqmcStats.h"
#include "src/core/io/ParallelStats.h"
#include "Propagator.h"

class Solver {

    size_t m_icycle = 0ul;
    Propagator &m_prop;
    const Options &m_opts;
    Wavefunction &m_wf;
    Reference m_reference;

    std::unique_ptr<FciqmcStats> m_stats = nullptr;
    std::unique_ptr<ParallelStats> m_parallel_stats = nullptr;

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

    conn::Basic<> m_connection;

    InteractiveVariable<bool> m_exit;
    std::unique_ptr<UniformTwf> m_uniform_twf;
    std::unique_ptr<StaticTwf> m_hubbard_twf;
    BilinearMevGroup m_mevs;

public:

    Solver(Propagator &prop, Wavefunction &wf, TableBase::Loc ref_loc);

    void execute(size_t niter=1);

    void begin_cycle();

    void propagate_row(const size_t& ipart);

    void loop_over_occupied_onvs();

    void annihilate_row(const size_t dst_ipart, const fields::Onv<>& dst_onv, const defs::wf_t& delta_weight, bool allow_initiation, const size_t& irow_store);

    void annihilate_row(const size_t dst_ipart, const fields::Onv<>& dst_onv, const defs::wf_t& delta_weight, bool allow_initiation) {
        annihilate_row(dst_ipart, dst_onv, delta_weight, allow_initiation, *m_wf.m_store[dst_onv]);
    }

    void make_mev_contribs(const fields::Onv<>& src_onv, const defs::wf_t& src_weight){
        // m_wf.m_store.m_row is assumed to have been moved to the store row of the dst ONV
        m_mevs.make_contribs(src_onv, src_weight, m_wf.m_store.m_row.m_onv, m_wf.m_store.m_row.m_weight[0]);
        //std::cout << src_onv.to_string() << " " << m_wf.m_store.m_row.m_onv.to_string() << std::endl;
    }

    void make_mev_contribs_from_unique_src_onvs(SpawnTableRow& row_current, SpawnTableRow& row_block_start,
                                                const size_t& irow_block_end, const size_t& irow_store){
        // if the dst onv is not stored, it cannot give contibutions to any MEVs
        if (irow_store==~0ul) {
            row_current.jump(irow_block_end);
            return;
        }
        m_wf.m_store.m_row.jump(irow_store);
        /*
         * similar approach to loop_over_spawned, except the "blocks" in this instance refer to groups
         * of contributions from the same source ONV. src_weights emitted by a stochastic propagator are
         * appropriately scaled by the probability that at least one excitation to dst_onv was drawn.
         */
        auto get_nrow_in_block = [&]() { return row_current.m_i - row_block_start.m_i; };
        row_block_start.jump(row_current);

        for (; row_current.m_i < irow_block_end; row_current.step()) {
            ASSERT(m_wf.m_store.m_row.m_onv == row_current.m_dst_onv);
            // seek to next "parent" ONV
            if (row_current.m_src_onv != row_block_start.m_src_onv) {
                ASSERT(get_nrow_in_block()>0);
                // row_current is pointing to the first row of the next src_onv block
                // row_block_start can be used to access the src ONV data
                make_mev_contribs(row_block_start.m_src_onv, row_block_start.m_src_weight);
                row_block_start.jump(row_current);
            }
        }
        // finish off last block
        make_mev_contribs(row_block_start.m_src_onv, row_block_start.m_src_weight);
    }

    void loop_over_spawned();

    void end_cycle();

    void output_stats();
};

#endif //M7_SOLVER_H
