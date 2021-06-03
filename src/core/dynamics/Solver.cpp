//
// Created by rja on 10/11/2020.
//

#include "Solver.h"

Solver::Solver(Propagator &prop, Wavefunction &wf, std::vector<TableBase::Loc> ref_locs) :
        m_prop(prop),
        m_opts(prop.m_opts),
        m_wf(wf),
        m_refs(m_opts, m_prop.m_ham, m_wf, ref_locs),
        m_connection(prop.m_ham.nsite()),
        m_exit("exit"),
        m_uniform_twf(m_opts.spf_uniform_twf ? new UniformTwf(m_wf.npart(), prop.m_ham.nsite()) : nullptr),
        m_weighted_twf(m_opts.spf_weighted_twf ?
                       new WeightedTwf(m_wf.npart(), prop.m_ham.nsite(),
                                       m_opts.spf_twf_fermion_factor,
                                       m_opts.spf_twf_boson_factor) :
                       nullptr),
        m_mevs(m_opts, prop.m_ham.nsite(), prop.m_ham.nelec()), m_detsub(m_wf) {

    if (defs::enable_mevs && m_opts.rdm_rank > 0) {

        if (!m_opts.replicate && !m_prop.is_exact() && !m_opts.mev_mixed_estimator)
            log::warn("Attempting a stochastic propagation estimation of MEVs without replication, "
                      "this is biased");

        if (m_opts.replicate && m_prop.is_exact())
            log::warn("Attempting an exact-propagation estimation of MEVs with replication, "
                                          "the replica population is redundant");

        if (!m_opts.replicate && !m_prop.is_exact() && !m_opts.mev_mixed_estimator)
            log::warn("Attempting a stochastic estimation of MEVs without replication, this is biased");

        if (m_opts.replicate && !m_opts.mev_mixed_estimator)
            log::warn("Attempting a mixed estimation of MEVs with replication, "
                      "the replica population is redundant");
    }

    if (mpi::i_am_root())
        m_stats = std::unique_ptr<FciqmcStats>(new FciqmcStats("M7.stats", "FCIQMC", {wf.m_format}));
    if (m_opts.parallel_stats)
        m_parallel_stats = std::unique_ptr<ParallelStats>(
                new ParallelStats("M7.stats." + std::to_string(mpi::irank()), "FCIQMC Parallelization", {}));
    m_wf.m_ra.activate(m_icycle);
}

void Solver::execute(size_t ncycle) {

    if (!m_opts.read_hdf5_fname.empty()) {
        hdf5::FileReader fr(m_opts.read_hdf5_fname);
        hdf5::GroupReader gr("solver", fr);
        //m_wf.h5_read(gr, m_prop.m_ham, m_reference.get_onv());
        //loop_over_spawned();
        if (m_mevs.m_fermion_rdm) {
            hdf5::GroupReader gr2("rdm", gr);
            m_mevs.m_fermion_rdm->h5_read(gr);
            m_mevs.m_fermion_rdm->end_cycle();
        }
    }

    for (size_t i = 0ul; i < ncycle; ++i) {
        m_cycle_timer.reset();
        m_cycle_timer.unpause();
        begin_cycle();

        m_propagate_timer.reset();
        m_propagate_timer.unpause();
        m_detsub.gather();
        loop_over_occupied_onvs();
        m_propagate_timer.pause();

        m_communicate_timer.reset();
        m_communicate_timer.unpause();
        m_wf.communicate();
        m_communicate_timer.pause();

        m_annihilate_timer.reset();
        m_annihilate_timer.unpause();
        loop_over_spawned();
        m_detsub.make_mev_contribs(m_mevs, m_refs[0].get_onv());
        m_detsub.project(m_prop.tau());
        m_annihilate_timer.pause();

        end_cycle();
        m_cycle_timer.pause();
        output_stats();
        output_mevs();
        ++m_icycle;

        if (m_exit.read() && m_exit.m_v) break;
        if (defs::enable_mevs && m_mevs.m_accum_epoch) {
            if (i == m_mevs.m_accum_epoch.icycle_start() + m_opts.ncycle_accumulate_mevs) break;
        }
    }

    if (defs::enable_mevs && m_mevs.m_accum_epoch) {
        // repeat the last cycle but do not perform any propagation
        --m_icycle;
        finalizing_loop_over_occupied_onvs();;
    }

    if (!m_opts.write_hdf5_fname.empty()) {
        hdf5::FileWriter fw(m_opts.write_hdf5_fname);
        //m_wf.h5_write(gw);
        hdf5::GroupWriter gw("solver", fw);
        if (m_mevs.m_fermion_rdm) {
            hdf5::GroupWriter gw2("rdm", gw);
            m_mevs.m_fermion_rdm->h5_write(gw2);
        }
    }
}

void Solver::begin_cycle() {
    m_chk_nwalker_local = m_wf.m_nwalker.m_local[{0, 0}] + m_wf.m_delta_nwalker.m_local[{0, 0}];
    m_chk_ninitiator_local = m_wf.m_ninitiator.m_local[{0, 0}] + m_wf.m_delta_ninitiator.m_local[{(0, 0)}];
    m_wf.begin_cycle();
    ASSERT(m_wf.m_nwalker.m_local[0] == 0);
    m_wf.m_ra.update(m_icycle);
    m_propagate_timer.reset();
    m_refs.begin_cycle();

    auto update_epoch = [&](const size_t &ncycle_wait) {
        if (m_opts.replicate) {
            if (m_prop.m_shift.m_variable_mode[0] && m_prop.m_shift.m_variable_mode[1]) {
                auto max_start = std::max(
                        m_prop.m_shift.m_variable_mode[0].icycle_start(),
                        m_prop.m_shift.m_variable_mode[1].icycle_start());
                if (m_icycle > max_start + ncycle_wait) return true;
            }
        } else {
            if (m_prop.m_shift.m_variable_mode[0]) {
                auto start = m_prop.m_shift.m_variable_mode[0].icycle_start();
                if (m_icycle > start + ncycle_wait) return true;
            }
        }
        return false;
    };

    if (defs::enable_mevs) {
        if (m_mevs.m_accum_epoch.update(m_icycle, update_epoch(m_opts.ncycle_wait_mevs))) {
            ASSERT(m_mevs.m_fermion_rdm->m_store.m_hwm == 0);
        }
    }

    if (m_opts.do_semistochastic && !m_detsub.m_epoch) {
        auto init = m_detsub.m_epoch.update(m_icycle, update_epoch(m_opts.ncycle_wait_detsub));
        if (init) {
            m_detsub.build_from_all_occupied(m_prop.m_ham);
            std::cout << m_wf.m_store.to_string() << std::endl;
            log::debug("initialized deterministic subspace");
        }
    }
}

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
    auto &row = m_wf.m_store.m_row;

    for (row.restart(); row.in_range(); row.step()) {
        /*
         * stats always refer to the state of the wavefunction in the previous iteration
         */

        if (row.m_onv.is_zero())
            /*
             * this is a free row, caused by an earlier call to m_wf.remove_row()
             */
            continue;


        if (row.m_weight.is_zero() && !row.is_protected()) {
            /*
             * ONV has become unoccupied in all parts and must be removed from mapped list, but it must first make all
             * associated averaged contributions to MEVs
             */
            make_average_weight_mev_contribs(m_icycle);
            m_wf.remove_row();
            continue;
        }

        /*
         * if the accumulation of MEVs has just started, treat the row as though it became occupied in the annihilation
         * loop of the last MC cycle.
         */
        if (defs::enable_mevs && m_mevs.m_accum_epoch.started_this_cycle(m_icycle)) {
            ASSERT(m_mevs.m_accum_epoch);
            row.m_icycle_occ = m_icycle-1;
            row.m_average_weight = 0;
        }

        if (defs::enable_mevs && m_mevs.m_accum_epoch) {
            row.m_average_weight += row.m_weight;
        }

        if (m_mevs.is_period_cycle(m_icycle)) {
            /*
             * this is the end of a planned block-averaging cycle, therefore there may be unaccounted-for contributions
             * which need to be included in the average
             */
            make_average_weight_mev_contribs(m_icycle);
        }

        for (size_t ipart = 0ul; ipart < m_wf.m_format.m_nelement; ++ipart) {

            MPI_ASSERT(!m_wf.m_store.m_row.m_onv.is_zero(),
                       "Stored ONV should not be zeroed");
            MPI_ASSERT(mpi::i_am(m_wf.get_rank(m_wf.m_store.m_row.m_onv)),
                       "Stored ONV should be on its allocated rank");

            const auto &weight = row.m_weight[ipart];

            m_wf.m_nocc_onv.m_local++;
            if (row.m_initiator.get(ipart))
                m_wf.m_ninitiator.m_local[ipart]++;

            m_wf.m_nwalker.m_local[ipart] += std::abs(weight);
            m_wf.m_l2_norm_square.m_local[ipart] += std::pow(std::abs(weight), 2.0);

            if (ipart == 0) {
                m_refs.contrib_row();
                if (m_opts.spf_uniform_twf) m_uniform_twf->add(m_prop.m_ham, row.m_weight, row.m_onv);
                if (m_opts.spf_weighted_twf) m_weighted_twf->add(m_prop.m_ham, row.m_weight, row.m_onv);
            }

            if (m_wf.m_ra.is_active()) {
                m_spawning_timer.reset();
                m_spawning_timer.unpause();
            }
            propagate_row(ipart);
            if (m_wf.m_ra.is_active()) {
                m_spawning_timer.pause();
                m_wf.m_ra.record_work_time(row, m_spawning_timer);
            }
        }
    }
    m_synchronization_timer.reset();
    m_synchronization_timer.unpause();
    mpi::barrier();
    m_synchronization_timer.pause();
}

void Solver::finalizing_loop_over_occupied_onvs() {
    if (!m_mevs.m_accum_epoch) return;
    auto &row = m_wf.m_store.m_row;
    for (row.restart(); row.in_range(); row.step()) {
        if (!row.m_onv.is_zero()) make_average_weight_mev_contribs(m_icycle);
    }
    if (m_mevs.m_fermion_rdm) m_mevs.m_fermion_rdm->end_cycle();
}

void Solver::annihilate_row(const size_t &dst_ipart, const fields::Onv<> &dst_onv, const defs::wf_t &delta_weight,
                            bool allow_initiation, bool src_deterministic, const size_t &irow_store) {
    if (m_opts.nadd_initiator == 0.0) ASSERT(allow_initiation);
    ASSERT(!dst_onv.is_zero());
    // check that the received determinant has come to the right place
    ASSERT(m_wf.m_ra.get_rank(dst_onv) == mpi::irank())
    // zero magnitude weights should not have been communicated
    if (consts::float_is_zero(delta_weight)) return;
    ASSERT(!consts::float_is_zero(delta_weight));

    m_wf.m_nspawned.m_local[dst_ipart] += delta_weight;
    if (irow_store == ~0ul) {
        /*
         * the destination ONV is not currently occupied, so initiator rules
         * must be applied
         */
        if (!allow_initiation) {
            //m_aborted_weight += std::abs(*delta_weight);
            return;
        }

        m_wf.create_row(m_icycle, dst_onv,m_prop.m_ham.get_energy(dst_onv), m_refs.is_connected(dst_onv));
        m_wf.set_weight(dst_ipart, delta_weight);

    } else {
        m_wf.m_store.m_row.jump(irow_store);
        if (src_deterministic && m_wf.m_store.m_row.m_deterministic.get(0)) return;
        defs::wf_t weight_before = m_wf.m_store.m_row.m_weight[dst_ipart];
        auto weight_after = weight_before + delta_weight;
        if (!consts::float_is_zero(weight_before) && !consts::float_is_zero(weight_after)
            && ((weight_before > 0) != (weight_after > 0)))
            m_wf.m_nannihilated.m_local[dst_ipart] += std::abs(std::abs(weight_before) - std::abs(weight_after));
        m_wf.change_weight(dst_ipart, delta_weight);
    }
}

void Solver::make_average_weight_mev_contribs(const size_t& icycle) {
    if (!defs::enable_mevs) return;
    if (!m_mevs.m_accum_epoch) return;
    auto &row = m_wf.m_store.m_row;
    // the current cycle should be included in the denominator
    auto ncycle_occ = row.occupied_ncycle(icycle);
    if (!ncycle_occ) {
        ASSERT(row.m_average_weight.is_zero());
        return;
    }

    for (size_t ipart = 0ul; ipart < m_wf.m_format.m_nelement; ++ipart) {
        auto& ref = m_refs[ipart];
        auto& ref_onv = ref.get_onv();
        auto ipart_replica = m_wf.ipart_replica(ipart);
        /*
         * if contributions are coming from two replicas, we should take the mean
         */
        double dupl_fac = (ipart_replica == ipart) ? 1.0 : 0.5;
        /*
         * the "average" weights actually refer to the unnormalized average. The averages are obtained by dividing
         * each by the number of cycles for which the row is occupied.
         */
        const auto av_weight = row.m_average_weight[ipart] / ncycle_occ;
        auto av_weight_rep = m_mevs.get_ket_weight(row.m_average_weight[ipart_replica] / ncycle_occ);
        /*
         * scale up the product by a factor of the number of instantaneous contributions being accounted for in this
         * single averaged contribution (ncycle_occ)
         */
        m_mevs.m_fermion_rdm->make_contribs(row.m_onv, dupl_fac*ncycle_occ*av_weight, row.m_onv, av_weight_rep);
        if (m_mevs.m_explicit_hf_conns) {
            if (row.m_reference_connection.get(0) && !ref.is_same(row)) {
                const auto av_weight_ref = ref.norm_average_weight(icycle, ipart);
                const auto av_weight_ref_rep = m_mevs.get_ket_weight(ref.norm_average_weight(icycle, ipart_replica));
                m_mevs.m_fermion_rdm->make_contribs(ref_onv, dupl_fac*ncycle_occ*av_weight_ref,
                                                    row.m_onv, av_weight_rep);
                m_mevs.m_fermion_rdm->make_contribs(row.m_onv, dupl_fac*ncycle_occ*av_weight,
                                                    ref_onv, av_weight_ref_rep);
            }
        }
    }
    row.m_average_weight = 0;
    row.m_icycle_occ = icycle;
}


void Solver::make_instant_mev_contribs(const fields::Onv<> &src_onv, const defs::wf_t &src_weight,
                                       const size_t &dst_ipart) {
    // m_wf.m_store.m_row is assumed to have jumped to the store row of the dst ONV
    if (!m_mevs.m_accum_epoch) return;
    if (m_mevs.m_explicit_hf_conns) {
        if (src_onv == m_refs[dst_ipart].get_onv() || m_wf.m_store.m_row.m_onv == m_refs[dst_ipart].get_onv()) return;
    }
    auto dst_ipart_replica = m_wf.ipart_replica(dst_ipart);
    double dupl_fac = (dst_ipart_replica == dst_ipart) ? 1.0 : 0.5;
    if (m_mevs.m_fermion_rdm) {
        /*
         * We need to be careful of the intermediate state of the walker weights.
         * if src_weight is taken from the wavefunction at cycle i, dst_weight is at an intermediate value equal to
         * the wavefunction at cycle i with the diagonal part of the propagator already applied. We don't want to
         * introduce a second post-annihilation loop over occupied ONVs to apply the diagonal part of the
         * propagator, so for MEVs, the solution is to reconstitute the value of the walker weight before the
         * diagonal death/cloning.
         *
         * The death-step behaviour of the exact (and stochastic on average) propagator is to scale the WF:
         * Ci -> Ci*(1 - tau (Hii-shift)).
         * By the time MEV contributions are being made, the death step has already been applied, and so the pre-
         * death value of the weight must be reconstituted by undoing the scaling, thus, the pre-death value of Cdst
         * is just Cdst/(1 - tau (Hii-shift))
         */
        auto dst_weight_before_death = m_wf.m_store.m_row.m_weight[dst_ipart_replica];
        dst_weight_before_death /=
                1 - m_prop.tau() * (m_wf.m_store.m_row.m_hdiag - m_prop.m_shift[dst_ipart_replica]);
        m_mevs.m_fermion_rdm->make_contribs(src_onv, dupl_fac * src_weight, m_wf.m_store.m_row.m_onv,
                                            dst_weight_before_death);
    }
}

void Solver::loop_over_spawned() {
    if (!m_wf.recv().m_hwm) return;
    mpi::barrier();
    if (m_opts.consolidate_spawns) {
        auto row1 = m_wf.recv().m_row;
        auto row2 = m_wf.recv().m_row;
        auto comp_fn = [&](const size_t &irow1, const size_t &irow2) {
            row1.jump(irow1);
            row2.jump(irow2);
            // sort criteria from major to minor: dst ONV, dst_ipart, src ONV,
            if (row1.m_dst_onv == row2.m_dst_onv) {
                if (row1.m_dst_ipart == row2.m_dst_ipart) {
                    return row1.m_src_onv <= row2.m_src_onv;
                }
                return row1.m_dst_ipart <= row2.m_dst_ipart;
            }
            return row1.m_dst_onv <= row2.m_dst_onv;
        };

        /*
         * sorting in ascending lexical order
         */
        QuickSorter qs(comp_fn);
        qs.reorder_sort(m_wf.recv());

        /*
         * now that the recv list is reordered according to dst ONV, we must process its
         * contents in dst ONV *blocks*.
         *
         * e.g.
         * irow  dst   src
         * 0     A     P
         * 1     A     P
         * 2     A     Q
         * 3     A     Q
         * 4     A     Q
         * 5     B     P
         * 6     B     R
         * 7     B     R
         * 8     C     Q
         * 9     C     Q
         * 10    C     R
         *
         * the algorithm uses two row indices:
         *      irow_block_start
         *      irow_current
         *
         * initially, both are 0.
         *
         * irow_current is iteratively incremented until its dst ONV does not equal that of
         * irow_block_start. at each iteration, the delta_weight is accumulated into a total.
         *
         * once this condition is met (at row 5 in the example), irow_block_start is iterated
         * to the row before irow_current,
         *
         * We also need to be careful of the walker weights
         * if Csrc is taken from the wavefunction at cycle i,
         * Cdst is at an intermediate value equal to the wavefunction at cycle i with the diagonal
         * part of the propagator already applied.
         *
         * We don't want to introduce a second post-annihilation loop over occupied ONVs to apply
         * the diagonal part of the propagator, so for MEVs, the solution is to reconstitute the
         * value of the walker weight before the diagonal death/cloning.
         *
         * In the exact propagator, the death step does:
         * Ci -> Ci*(1 - tau (Hii-shift)).
         *
         * thus, the pre-death value of Cdst is just Cdst/(1 - tau (Hii-shift))
         *
         */

        auto row_block_start = m_wf.recv().m_row;
        auto row_current = m_wf.recv().m_row;
        auto row_block_start_src_blocks = m_wf.recv().m_row;

        auto get_nrow_in_block = [&]() { return row_current.m_i - row_block_start.m_i; };
        auto get_allow_initiation = [&]() {
            // row_block_start is now at last row in last block
            bool allow = get_nrow_in_block() > 1;
            if (!allow) {
                // only one src_onv for this dst_onv. If the parent is an initiator,
                // contributions to unoccupied ONVs are allowed
                allow = row_block_start.m_src_initiator;
            }
            return allow;
        };

        auto still_in_block = [&]() {
            return row_current.m_dst_onv == row_block_start.m_dst_onv &&
                   row_current.m_dst_ipart == row_block_start.m_dst_ipart;
        };

        row_block_start.restart();
        defs::wf_t total_delta = 0.0;
        for (row_current.restart(); row_current.in_range(); row_current.step()) {
            size_t dst_ipart = row_block_start.m_dst_ipart;
            if (still_in_block()) {
                total_delta += row_current.m_delta_weight;
            } else {
                // row_current is in first row of next block
                ASSERT(get_nrow_in_block() > 0);
                // get the row index (if any) of the dst_onv
                auto irow_store = *m_wf.m_store[row_block_start.m_dst_onv];
                make_mev_contribs_from_unique_src_onvs(row_block_start, row_block_start_src_blocks,
                                                       row_current.m_i - 1, irow_store);
                ASSERT(row_block_start.m_i == row_current.m_i - 1)
                annihilate_row(dst_ipart, row_block_start.m_dst_onv, total_delta, get_allow_initiation(), irow_store);
                // put block start to start of next block
                row_block_start.step();
                ASSERT(row_block_start.m_i == row_current.m_i)
                total_delta = row_current.m_delta_weight;
            }
        }
        // finish off last block
        if (row_block_start.in_range()) {
            size_t dst_ipart = row_block_start.m_dst_ipart;
            auto irow_store = *m_wf.m_store[row_block_start.m_dst_onv];
            make_mev_contribs_from_unique_src_onvs(row_block_start, row_block_start_src_blocks,
                                                   m_wf.recv().m_hwm - 1, irow_store);
            annihilate_row(dst_ipart, row_block_start.m_dst_onv, total_delta,
                           get_allow_initiation(), row_block_start.m_src_deterministic, irow_store);
        }
    } else {
        auto &row = m_wf.recv().m_row;
        if (m_opts.rdm_rank > 0 && m_mevs.m_accum_epoch) {
            /*
             * an additional loop over recvd spawns is required in this case, in order to make the necessary
             * MEV contributions before the new spawns are added to the instantaneous populations
             */
            for (row.restart(); row.in_range(); row.step()) {
                auto irow_store = *m_wf.m_store[row.m_dst_onv];
                if (irow_store != ~0ul) {
                    m_wf.m_store.m_row.jump(irow_store);
                    if (row.m_src_deterministic && m_wf.m_store.m_row.m_deterministic.get(0)) {
                        continue;
                    }
                    make_instant_mev_contribs(row.m_src_onv, row.m_src_weight, row.m_dst_ipart);
                }
            }
        }

        for (row.restart(); row.in_range(); row.step()) {
            annihilate_row(row.m_dst_ipart, row.m_dst_onv, row.m_delta_weight, row.m_src_initiator,
                           row.m_src_deterministic);
        }
    }
    m_wf.recv().clear();
}

void Solver::propagate_row(const size_t &ipart) {
    auto &row = m_wf.m_store.m_row;

    if (row.is_cleared()) return;

    if (consts::float_is_zero(row.m_weight[ipart])) return;

    m_prop.off_diagonal(m_wf, ipart);
    m_prop.diagonal(m_wf, ipart);
}

void Solver::end_cycle() {
    /*
     * TODO: make these checks compatible with dynamic rank allocation
     */
//    double chk_ratio;
    if (!consts::float_is_zero(m_wf.m_nwalker.m_local[{0, 0}])) {
//        chk_ratio = m_chk_nwalker_local / m_wf.m_nwalker(0, 0);
//        bool chk = m_chk_nwalker_local == 0.0 || consts::floats_nearly_equal(chk_ratio, 1.0);
//        if (!chk) std::cout << "discrepancy: " << m_chk_nwalker_local-m_wf.m_nwalker(0, 0) << std::endl;
//        MPI_REQUIRE(chk,"Unlogged walker population changes have occurred");
    }
//    MPI_REQUIRE(m_chk_ninitiator_local == m_wf.m_ninitiator(0, 0),
//                "Unlogged creations of initiator ONVs have occurred");

    if (m_mevs.m_fermion_rdm) m_mevs.m_fermion_rdm->end_cycle();
    m_wf.end_cycle();
    m_refs.end_cycle();
    m_prop.update(m_icycle, m_wf);
    if (m_uniform_twf) m_uniform_twf->reduce();
    if (m_weighted_twf) m_weighted_twf->reduce();
}

void Solver::output_stats() {

    auto sync_overhead = mpi::all_sum((double) m_synchronization_timer);
    if (mpi::i_am_root()) {
        auto &stats = m_stats->m_row;
        stats.m_icycle = m_icycle;
        stats.m_tau = m_prop.tau();
        stats.m_shift = m_prop.m_shift.m_values;
        stats.m_nwalker = m_wf.m_nwalker.m_reduced;
        stats.m_delta_nwalker = m_wf.m_delta_nwalker.m_reduced;
        stats.m_nwalker_spawned = m_wf.m_nspawned.m_reduced;
        stats.m_nwalker_annihilated = m_wf.m_nannihilated.m_reduced;
        stats.m_ref_proj_energy_num = m_refs[0].proj_energy_num();
        stats.m_ref_weight = m_refs[0].weight();
        stats.m_ref_proj_energy = stats.m_ref_proj_energy_num;
        stats.m_ref_proj_energy /= stats.m_ref_weight;
        stats.m_l2_norm = m_wf.m_l2_norm_square.m_reduced;
        stats.m_l2_norm.to_sqrt();
        stats.m_ninitiator = m_wf.m_ninitiator.m_reduced;
        stats.m_nocc_onv = m_wf.m_nocc_onv.m_reduced;
        stats.m_delta_nocc_onv = m_wf.m_delta_nocc_onv.m_reduced;
        stats.m_psingle = m_prop.m_magnitude_logger.m_psingle;
        stats.m_total_synchronization_overhead = sync_overhead;
        stats.m_propagate_loop_time = m_propagate_timer;
        stats.m_communication_time = m_communicate_timer;
        stats.m_annihilation_loop_time = m_annihilate_timer;
        stats.m_total_cycle_time = m_cycle_timer;
        if (m_uniform_twf) stats.m_uniform_twf_num = m_uniform_twf->m_numerator_total[0];
        if (m_weighted_twf) {
            stats.m_weighted_twf_num = m_weighted_twf->m_numerator_total[0];
            stats.m_weighted_twf_denom = m_weighted_twf->m_denominator_total[0];
        }
        stats.m_reweighting_factor = m_prop.m_shift.m_reweighter.m_total;
        m_stats->flush();
    }

    if (m_opts.parallel_stats) {
        auto &stats = m_parallel_stats->m_row;
        stats.m_icycle = m_icycle;
        stats.m_synchronization_overhead = m_synchronization_timer;
        stats.m_nblock_wf_ra = m_wf.m_ra.nblock_local();
        stats.m_nwalker_total = m_wf.m_nwalker.m_reduced.sum();
        stats.m_nwalker_lookup_skip = m_wf.m_store.m_ntotal_skip;
        stats.m_nwalker_lookup = m_wf.m_store.m_ntotal_lookup;
        stats.m_nrow_recv = m_wf.m_comm.m_last_recv_count;
        m_parallel_stats->flush();
    }
}

void Solver::annihilate_row(const size_t &dst_ipart, const fields::Onv<> &dst_onv, const defs::wf_t &delta_weight,
                            bool allow_initiation, bool src_deterministic) {
    annihilate_row(dst_ipart, dst_onv, delta_weight, allow_initiation, src_deterministic, *m_wf.m_store[dst_onv]);
}

void Solver::make_mev_contribs_from_unique_src_onvs(SpawnTableRow &row_current, SpawnTableRow &row_block_start,
                                                    const size_t &irow_block_end, const size_t &irow_store) {
    // if the dst onv is not stored, it cannot give contributions to any MEVs
    if (irow_store == ~0ul) {
        row_current.jump(irow_block_end);
        return;
    }
    m_wf.m_store.m_row.jump(irow_store);
    /*
     * similar approach to loop_over_spawned, except the "blocks" in this instance refer to groups
     * of contributions from the same source ONV. src_weights emitted by a stochastic propagator are
     * appropriately scaled by the probability that at least one excitation to dst_onv was drawn.
     */

    row_block_start.jump(row_current);

    for (; row_current.in_range(irow_block_end); row_current.step()) {
        ASSERT(m_wf.m_store.m_row.m_onv == row_current.m_dst_onv);
        // seek to next "parent" ONV
        if (row_current.m_src_onv != row_block_start.m_src_onv) {
            ASSERT(row_current.m_i - row_block_start.m_i > 0);
            // row_current is pointing to the first row of the next src_onv block
            // row_block_start can be used to access the src ONV data
            make_instant_mev_contribs(row_block_start.m_src_onv, row_block_start.m_src_weight,
                                      row_block_start.m_dst_ipart);
            row_block_start.jump(row_current);
        }
    }
    // finish off last block
    make_instant_mev_contribs(row_block_start.m_src_onv, row_block_start.m_src_weight, row_block_start.m_dst_ipart);
}

void Solver::output_mevs() {
    if (!(defs::enable_mevs && m_mevs.is_period_cycle(m_icycle))) return;
    hdf5::FileWriter fw(std::to_string(m_mevs.iperiod(m_icycle)) + "." + m_opts.write_hdf5_fname);
    hdf5::GroupWriter gw("solver", fw);
    if (m_mevs.m_fermion_rdm) {
        hdf5::GroupWriter gw2("rdm", gw);
        m_mevs.m_fermion_rdm->h5_write(gw2);
        const auto b = -99.9421389039331;
        std::cout << m_mevs.m_fermion_rdm->get_energy(m_prop.m_ham)-b << " "
                  << m_prop.m_shift.m_values[0]-b << " " << m_refs[0].proj_energy_num()[0]/m_refs[0].weight()[0]-b << std::endl;
    }
}

const MevGroup &Solver::mevs() const {
    return m_mevs;
}
