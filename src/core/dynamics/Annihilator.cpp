//
// Created by rja on 14/07/2021.
//

#include "Annihilator.h"

comparators::index_cmp_fn_t Annihilator::make_sort_cmp_fn() {
    if (m_rdms) {
        return [&](const size_t &irow1, const size_t &irow2) {
            m_work_row1.jump(irow1);
            m_work_row2.jump(irow2);
            // sort criteria from major to minor: dst ONV, dst_ipart, src ONV,
            if (m_work_row1.m_dst_mbf == m_work_row2.m_dst_mbf) {
                if (m_work_row1.m_ipart_dst == m_work_row2.m_ipart_dst) {
                    return m_work_row1.m_src_mbf < m_work_row2.m_src_mbf;
                }
                return m_work_row1.m_ipart_dst < m_work_row2.m_ipart_dst;
            }
            return m_work_row1.m_dst_mbf < m_work_row2.m_dst_mbf;
        };
    } else {
        return [&](const size_t &irow1, const size_t &irow2) {
            m_work_row1.jump(irow1);
            m_work_row2.jump(irow2);
            // sort criteria from major to minor: dst ONV, dst_ipart, src ONV,
            if (m_work_row1.m_dst_mbf == m_work_row2.m_dst_mbf) {
                return m_work_row1.m_ipart_dst < m_work_row2.m_ipart_dst;
            }
            return m_work_row1.m_dst_mbf < m_work_row2.m_dst_mbf;
        };
    }
}

Annihilator::Annihilator(Wavefunction &wf, const Propagator &prop, const References &refs,
                         Rdms &rdms, const size_t &icycle, defs::wf_comp_t nadd) :
        m_wf(wf), m_prop(prop), m_refs(refs), m_rdms(rdms), m_nadd(nadd), m_icycle(icycle),
        m_work_row1(wf.m_comm.recv().m_row), m_work_row2(wf.m_comm.recv().m_row),
        m_sort_cmp_fn(make_sort_cmp_fn()) {
    REQUIRE_TRUE_ALL(bool(m_rdms)==wf.m_comm.recv().m_row.m_send_parents,
                     "cannot sample RDMs through annihilation unless parent MBFs are communicated");
}

void Annihilator::sort_recv() {
    QuickSorter qs(m_sort_cmp_fn);
    qs.reorder_sort(m_wf.recv());
}

void Annihilator::annihilate_row(const size_t &dst_ipart, const field::Mbf &dst_mbf, const defs::wf_t &delta_weight,
                                 bool allow_initiation, WalkerTableRow &dst_row) {
    if (m_nadd == 0.0) {
        DEBUG_ASSERT_TRUE(allow_initiation,
                          "initiator rules are turned off, every initiating annihilation should be allowed");
    }
    DEBUG_ASSERT_FALSE(dst_mbf.is_zero(), "recv table is contiguous, shouldn't have any cleared rows");
    DEBUG_ASSERT_EQ(m_wf.m_ra.get_rank(dst_mbf), mpi::irank(),
                    "the received ONV has not been sent to the right MPI rank!")
    // zero magnitude weights should not have been communicated
    if (consts::float_is_zero(delta_weight)) return;

    m_wf.m_nspawned.m_local[dst_ipart] += std::abs(delta_weight);
    if (!dst_row.in_range()) {
        /*
         * the destination MBF row in m_wf.m_store is not currently occupied, so initiator rules
         * must be applied
         */
        if (!allow_initiation) {
            //m_aborted_weight += std::abs(*delta_weight);
            return;
        }

        m_wf.create_row_(m_icycle, dst_mbf, m_prop.m_ham.get_energy(dst_mbf), m_refs.is_connected(dst_mbf));
        m_wf.set_weight(dst_ipart, delta_weight);

    } else {
        defs::wf_t weight_before = dst_row.m_weight[dst_ipart];
        if (consts::float_is_zero(weight_before) && !allow_initiation) {
            // the row exists, but that does not necessarily mean that each part has non-zero occupation
            //m_aborted_weight += std::abs(*delta_weight);
            return;
        }
        auto weight_after = weight_before + delta_weight;
        if (!consts::float_is_zero(weight_before) && !consts::float_is_zero(weight_after)
            && ((std::abs(weight_before) > 0) != (std::abs(weight_after) > 0)))
            m_wf.m_nannihilated.m_local[dst_ipart] += std::abs(std::abs(weight_before) - std::abs(weight_after));
        m_wf.change_weight(dst_ipart, delta_weight);
    }
}

void Annihilator::handle_dst_block(SpawnTableRow &block_begin, SpawnTableRow &next_block_begin,
                                   const defs::wf_t &total_delta, WalkerTableRow &dst_row) {
    DEBUG_ASSERT_FALSE(in_same_dst_block(block_begin, next_block_begin),
                       "start of block and start of next block should not be in the same block");
    DEBUG_ASSERT_LT(block_begin.index(), next_block_begin.index(),
                    "block start should be strictly before current row in recv table");

    /*
     * only consider RDM contributions if there are any RDMs being accumulated and the destination exists
     */
    if (m_rdms.m_accum_epoch && dst_row.in_range()) {
        /*
         * store the original positions of the row objects in the recv table
         */
        size_t irow_begin = block_begin.index();
        size_t irow_next_begin = next_block_begin.index();

        auto &current = next_block_begin;
        defs::wf_t src_weight = block_begin.m_src_weight;
        DEBUG_ONLY(src_weight);

        for (current.jump(block_begin);; current.step()) {
            if (!in_same_src_block(current, block_begin)) {
                handle_src_block(block_begin, dst_row);
                block_begin.jump(current);
                /*
                 * the above call iterates block_begin through the contributions until the row index matches that of current
                 * which would be the start of the next block if current is not past the end of recv rows. In any case, the
                 * two rows should be the same. If the new start of block is past the end of recv rows, we're done
                 */
                if (block_begin.index() == irow_next_begin) break;
                DEBUG_ASSERT_EQ(block_begin.index(), current.index(),
                                "block_begin should have been pointed to the beginning of the next src_mbf block");
                src_weight = block_begin.m_src_weight;
            } else {
                DEBUG_ASSERT_EQ(current.m_dst_mbf, block_begin.m_dst_mbf, "dst MBFs should be the same");
                DEBUG_ASSERT_EQ(current.m_ipart_dst, block_begin.m_ipart_dst, "dst iparts should be the same");
                DEBUG_ASSERT_EQ(current.m_src_mbf, block_begin.m_src_mbf, "src MBFs should be the same");
                DEBUG_ASSERT_TRUE(consts::floats_equal(src_weight, defs::wf_t(current.m_src_weight)),
                                  "all spawns with the same dst_mbf, ipart_dst, and src_mbf should have the same src_weight");
            }
        }
        /*
         * restore the original positions
         */
        block_begin.jump(irow_begin);
        next_block_begin.jump(irow_next_begin);
    }

    bool allow_initiation = (next_block_begin.index() - block_begin.index()) > 1;
    if (!allow_initiation) {
        // only one src_mbf for this dst_mbf. If the parent is an initiator,
        // contributions to unoccupied ONVs are allowed
        allow_initiation = block_begin.m_src_initiator;
    }
    annihilate_row(0, block_begin.m_dst_mbf, total_delta, allow_initiation, dst_row);
    block_begin.jump(next_block_begin);
}

void Annihilator::handle_src_block(SpawnTableRow &block_begin, WalkerTableRow &dst_row) {
    DEBUG_ASSERT_TRUE(m_rdms.m_accum_epoch, "shouldn't be sampling RDMs yet");
    DEBUG_ASSERT_EQ(block_begin.m_dst_mbf, dst_row.m_mbf, "wrong dst_row found");

    size_t ipart_dst = block_begin.m_ipart_dst;

    /*
     * don't make contributions to RDM elements if they already take the equivalent contribution from deterministic
     * average connections to the reference
     */
    if (m_rdms.m_explicit_ref_conns){
        if (dst_row.m_mbf == m_refs[ipart_dst].get_mbf()) return;
        if (block_begin.m_src_mbf == m_refs[m_wf.ipart_replica(ipart_dst)].get_mbf()) return;
    }
    /*
     * or if they already take contributions from deterministic subspace connections.
     */
    if (block_begin.m_src_deterministic && dst_row.m_deterministic.get(ipart_dst)) return;

    m_rdms.make_contribs(block_begin, dst_row, m_prop);
}

void Annihilator::loop_over_dst_mbfs() {
    if (!m_wf.recv().m_hwm) return;

    /*
     * put the block_begin row to the beginning of the recv table
     */
    auto &block_begin = m_work_row1;
    block_begin.restart();

    /*
     * we need to be able to lookup the dst in the main table. these variables are for
     */
    DstFinder dst_finder(m_wf, block_begin);

    defs::wf_t total_delta = 0.0;
    dst_finder();

    auto &current = m_work_row2;
    for (current.restart();; current.step()) {
        if (!in_same_dst_block(current, block_begin)) {
            /*
             * either the current row has been iterated over the high water mark of the recv table, or onto a row with a
             * different (dst_mbf, ipart_dst) pair. In either case, we have reached the end of a block of the major
             * sorting field, so we must handle the block just finished
             */
            handle_dst_block(block_begin, current, total_delta, dst_finder.store_row());
            /*
             * the above call iterates block_begin through the contributions until the row index matches that of current
             * which would be the start of the next block if current is not past the end of recv rows. In any case, the
             * two rows should be the same. If the new start of block is past the end of recv rows, we're done
             */
            if (!block_begin.in_range()) break;
            DEBUG_ASSERT_EQ(block_begin.index(), current.index(),
                            "block_begin should have been pointed to the beginning of the next block");
            /*
             * new block of dst_mbfs, so reset total contrib
             */
            total_delta = 0.0;
            dst_finder();
        } else {
            DEBUG_ASSERT_NE(current.m_dst_mbf, current.m_src_mbf,
                            "should never have diagonal connections at annihilation");
            DEBUG_ASSERT_EQ(current.m_dst_mbf, block_begin.m_dst_mbf, "dst MBFs should be the same");
            DEBUG_ASSERT_EQ(current.m_ipart_dst, block_begin.m_ipart_dst, "dst iparts should be the same");
        }
        if (!dst_finder.m_deterministic || !current.m_src_deterministic) {
            // this is not a determ->determ connection, so include it
            total_delta += current.m_delta_weight;
        }
    }
}
