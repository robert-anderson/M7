//
// Created by rja on 10/11/2020.
//

#include "Solver.h"

void Solver::loop_over_occupied_onvs() {
    /*
     * Loop over all rows in the m_wf.m_walkers table and if the row is not empty:
     *      ascertain the initiator status of the ONV
     *      if not initiator and weight > initiator threshold, then grant initiator status
     *      if initiator and weight < initiator threshold, then revoke initiator status
     *      perform all the off-diagonal propagation (fill send table)
     *      update local weight in the diagonal cloning/death step
     * else if all elements of the m_weight field are zero, the row should be removed
     */
    auto& row = m_wf.m_store.m_row;
    for (row.restart(); row.in_range(); row.step()){
        /*
         * stats always refer to the state of the wavefunction in the previous iteration
         */

        if (row.is_cleared()) continue;

        MPI_ASSERT(!m_wf.m_store.m_row.m_onv.is_zero(),
                   "Stored ONV should not be zeroed");
        MPI_ASSERT(mpi::i_am(m_wf.get_rank(m_wf.m_store.m_row.m_onv)),
                   "Stored ONV should be on its allocated rank");

        const auto &weight = row.m_weight(m_wf.m_ipart);

        m_wf.m_nocc_onv(0, 0)++;
        if (row.m_initiator.get(m_wf.m_ipart))
            m_wf.m_ninitiator(0, 0)++;

        if (consts::float_is_zero(weight)){
            // ONV has become unoccupied and must be removed from mapped list
            m_wf.remove_walker();
            continue;
        }

        m_wf.m_nwalker(0, 0) += std::abs(weight);
        m_wf.m_l2_norm_square(0, 0) += std::pow(std::abs(weight), 2.0);

        m_reference.add_row();

        /*
        if (m_prop.m_variable_shift){
            m_connection.connect(m_reference.get_onv(), row.m_onv);
            if (m_connection.nexcit()==0) m_average_coeffs.m_ref_coeff += weight;
            if (m_connection.nexcit()==2) {
                auto irow = *m_average_coeffs[m_connection];
                if (irow==~0ul) irow = m_average_coeffs.insert(m_connection);
                m_average_coeffs.m_row.jump(irow);
                m_average_coeffs.m_row.m_values(0)+=weight;
            }
        }
         */

        if (m_wf.m_ra.is_active()) {
            m_spawning_timer.reset();
            m_spawning_timer.unpause();
        }
        propagate_row();
        if (m_wf.m_ra.is_active()) {
            m_spawning_timer.pause();
            m_wf.m_ra.record_work_time(row, m_spawning_timer);
        }
    }
    m_synchronization_timer.reset();
    m_synchronization_timer.unpause();
    mpi::barrier();
    m_synchronization_timer.pause();
}

void Solver::annihilate_row() {
    const auto& row = m_wf.recv().m_row;
    const auto& dst_onv = row.m_dst_onv;
    ASSERT(!dst_onv.is_zero());
    // check that the received determinant has come to the right place
    ASSERT(m_wf.m_ra.get_rank(dst_onv) == mpi::irank())
    const defs::wf_t& delta_weight = row.m_delta_weight;
    // zero magnitude weights should not have been communicated
    ASSERT(!consts::float_is_zero(delta_weight));

    auto irow_walkers = *m_wf.m_store[dst_onv];

    if (irow_walkers == ~0ul) {
        /*
         * the destination ONV is not currently occupied, so initiator rules
         * must be applied
         */
        if (!row.m_src_initiator) {
            //m_aborted_weight += std::abs(*delta_weight);
            return;
        }

        m_wf.create_walker_(
                dst_onv,
                delta_weight,
                m_prop.m_ham.get_energy(dst_onv),
                m_reference.is_connected(dst_onv));
    } else {
        m_wf.m_store.m_row.jump(irow_walkers);
        defs::wf_t weight_before = m_wf.m_store.m_row.m_weight(0);
        auto weight_after = weight_before+delta_weight;
        if ((weight_before>0)!=(weight_after>0))
            m_wf.m_nannihilated(0, 0) += std::abs(std::abs(weight_before)-std::abs(weight_after));
        m_wf.change_weight(delta_weight);
    }
}

void Solver::loop_over_spawned() {
    mpi::barrier();
    const auto &row = m_wf.recv().m_row;
    for (row.restart(); row.in_range(); row.step()){
        annihilate_row();
    }
    m_wf.recv().clear();
}

Solver::Solver(Propagator &prop, Wavefunction &wf, TableBase::Loc ref_loc) :
        m_prop(prop),
        m_opts(prop.m_opts),
        m_wf(wf),
        m_reference(m_opts, m_prop.m_ham, m_wf, 0, ref_loc),
        m_connection(prop.m_ham.nsite())
        //m_average_coeffs("average coeffs", {2, 2}, 1)
        {
    if (mpi::i_am_root())
        m_stats = mem_utils::make_unique<StatsFile<FciqmcStatsSpecifier>>("M7.stats");
    if (m_opts.parallel_stats)
        m_parallel_stats = mem_utils::make_unique<StatsFile<ParallelStatsSpecifier>>(
            "M7.stats."+std::to_string(mpi::irank()));
    m_wf.m_ra.activate(m_icycle);
}

void Solver::execute(size_t niter) {
    for (size_t i = 0ul; i < niter; ++i) {

        m_cycle_timer.reset();
        m_cycle_timer.unpause();
        begin_cycle();

        m_propagate_timer.reset();
        m_propagate_timer.unpause();
        loop_over_occupied_onvs();
        m_propagate_timer.pause();

        m_communicate_timer.reset();
        m_communicate_timer.unpause();
        m_wf.communicate();
        m_communicate_timer.pause();

        m_annihilate_timer.reset();
        m_annihilate_timer.unpause();
        //consolidate_spawned();
        loop_over_spawned();
        m_annihilate_timer.pause();

        end_cycle();
        m_cycle_timer.pause();
        output_stats();
        ++m_icycle;
    }
}

void Solver::propagate_row() {
    auto& row = m_wf.m_store.m_row;
    const auto& ipart = m_wf.m_ipart;

    if (row.is_cleared()) return;

    if (consts::float_is_zero(row.m_weight(ipart))) return;

    //bool is_deterministic = row.m_deterministic.get(ipart);

    m_prop.off_diagonal(m_wf);
    m_prop.diagonal(m_wf);
}

void Solver::begin_cycle() {
    m_chk_nwalker_local = m_wf.m_nwalker(0, 0) + m_wf.m_delta_nwalker(0, 0);
    m_chk_ninitiator_local = m_wf.m_ninitiator(0, 0) + m_wf.m_delta_ninitiator(0, 0);
    m_wf.begin_cycle();
    ASSERT(m_wf.m_nwalker(0,0)==0);
    m_wf.m_ra.update(m_icycle);
    m_propagate_timer.reset();
    m_reference.begin_cycle();
}

void Solver::end_cycle() {
    /*
     * TODO: make these checks compatible with dynamic rank allocation
     */
//    double chk_ratio;
    if (!consts::float_is_zero(m_wf.m_nwalker(0, 0))) {
//        chk_ratio = m_chk_nwalker_local / m_wf.m_nwalker(0, 0);
//        bool chk = m_chk_nwalker_local == 0.0 || consts::floats_nearly_equal(chk_ratio, 1.0);
//        if (!chk) std::cout << "discrepancy: " << m_chk_nwalker_local-m_wf.m_nwalker(0, 0) << std::endl;
//        MPI_REQUIRE(chk,"Unlogged walker population changes have occurred");
    }
//    MPI_REQUIRE(m_chk_ninitiator_local == m_wf.m_ninitiator(0, 0),
//                "Unlogged creations of initiator ONVs have occurred");

    m_wf.end_cycle();
    m_reference.end_cycle();
    m_prop.update(m_icycle, m_wf);
}

void Solver::output_stats() {

    auto sync_overhead = mpi::all_sum((double)m_synchronization_timer);
    if (mpi::i_am_root()) {
        m_stats->m_icycle() = m_icycle;
        m_stats->m_tau() = m_prop.tau();
        m_stats->m_shift() = m_prop.m_shift;
        m_stats->m_nwalker() = m_wf.m_nwalker.reduced(0, 0);
        m_stats->m_delta_nwalker() = m_wf.m_delta_nwalker.reduced(0, 0);
        m_stats->m_nwalker_annihilated() = m_wf.m_nannihilated.reduced(0, 0);
        m_stats->m_ref_proj_energy_num() = m_reference.proj_energy_num();
        m_stats->m_ref_weight() = m_reference.get_weight();
        m_stats->m_ref_proj_energy() = m_reference.proj_energy();
        m_stats->m_l2_norm() = std::sqrt(m_wf.m_l2_norm_square.reduced(0, 0));
        m_stats->m_ninitiator() = m_wf.m_ninitiator.reduced(0, 0);
        m_stats->m_nocc_onv() = m_wf.m_nocc_onv.reduced(0, 0);
        m_stats->m_delta_nocc_onv() = m_wf.m_delta_nocc_onv.reduced(0, 0);
        m_stats->m_psingle() = m_prop.m_magnitude_logger.m_psingle;
        m_stats->m_total_synchronization_overhead() = sync_overhead;
        m_stats->m_propagate_loop_time() = m_propagate_timer;
        m_stats->m_communication_time() = m_communicate_timer;
        m_stats->m_annihilation_loop_time() = m_annihilate_timer;
        m_stats->m_total_cycle_time() = m_cycle_timer;
        m_stats->flush();
    }

    if (m_opts.parallel_stats){
        m_parallel_stats->m_icycle() = m_icycle;
        m_parallel_stats->m_synchronization_overhead() = m_synchronization_timer;
        m_parallel_stats->m_nblock_wf_ra() = m_wf.m_ra.nblock_();
        m_parallel_stats->m_nwalker() = m_wf.m_nwalker(0, 0);
        m_parallel_stats->m_nwalker_lookup_skip() = m_wf.m_store.m_ntotal_skip;
        m_parallel_stats->m_nwalker_lookup() = m_wf.m_store.m_ntotal_lookup;
//    m_parallel_stats->m_nrow_free_walker_list() = m_wf.m_walkers.
//    StatsColumn<size_t> m_walker_list_high_water_mark;
//    StatsColumn<double> m_walker_list_high_water_mark_fraction;
//    StatsColumn<size_t> m_nrow_sent;
//    StatsColumn<size_t> m_largest_nrow_sent;
//    StatsColumn<double> m_largest_send_list_filled_fraction;
//    StatsColumn<size_t> m_irank_largest_nrow_sent;
    m_parallel_stats->m_nrow_recv() = m_wf.m_comm.m_last_recv_count;
//    StatsColumn<double> m_recv_list_filled_fraction;
        m_parallel_stats->flush();
    }
}
