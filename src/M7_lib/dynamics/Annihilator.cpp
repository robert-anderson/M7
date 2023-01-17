//
// Created by Robert J. Anderson on 14/07/2021.
//

#include "Annihilator.h"

comparators::index_cmp_fn_t Annihilator::make_sort_cmp_fn() {
    if (m_rdms) {
        return [&](const uint_t &irow1, const uint_t &irow2) {
            m_work_row1.jump(irow1);
            m_work_row2.jump(irow2);
            // sort criteria from major to minor: dst MBF, dst_ipart, src MBF,
            if (m_work_row1.m_dst_mbf == m_work_row2.m_dst_mbf) {
                if (m_work_row1.m_ipart_dst == m_work_row2.m_ipart_dst) {
                    return m_work_row1.m_src_mbf < m_work_row2.m_src_mbf;
                }
                return m_work_row1.m_ipart_dst < m_work_row2.m_ipart_dst;
            }
            return m_work_row1.m_dst_mbf < m_work_row2.m_dst_mbf;
        };
    } else {
        return [&](const uint_t &irow1, const uint_t &irow2) {
            m_work_row1.jump(irow1);
            m_work_row2.jump(irow2);
            // sort criteria from major to minor: dst MBF, dst_ipart
            if (m_work_row1.m_dst_mbf == m_work_row2.m_dst_mbf) {
                return m_work_row1.m_ipart_dst < m_work_row2.m_ipart_dst;
            }
            return m_work_row1.m_dst_mbf < m_work_row2.m_dst_mbf;
        };
    }
}

Annihilator::Annihilator(Wavefunction &wf, const Propagator &prop, const References &refs,
                         const shared_rows::Walker* hf, Rdms &rdms, const uint_t &icycle, wf_comp_t nadd) :
        m_wf(wf), m_prop(prop), m_refs(refs), m_hf(hf), m_rdms(rdms), m_nadd(nadd), m_icycle(icycle),
        m_work_row1(wf.m_send_recv.m_row), m_work_row2(wf.m_send_recv.m_row), m_dst_walker(m_wf.m_store.m_row),
        m_dst_weight(m_wf.npart(), wf_t{}), m_sort_cmp_fn(make_sort_cmp_fn()) {
    REQUIRE_TRUE_ALL(bool(m_rdms)==m_work_row1.m_send_parents,
                     "cannot sample RDMs through annihilation unless parent MBFs are communicated");
}

void Annihilator::sort_recv() {
    quicksort::Sorter qs(m_sort_cmp_fn);
    qs.reorder_sort(m_wf.recv());
}

void Annihilator::lookup_dst(const Mbf &dst_mbf, uint_t ipart_dst, bool &deterministic) {
    /*
     * if the current walker is valid and has MBF equal to dst_mbf, then this is not the first part being
     * contributed to in this annihilation loop, and so m_dst_walker and m_dst_weight should be left as-is
     */
    const auto same_dst_again = m_dst_walker && (m_dst_walker.m_mbf == dst_mbf);
    if (!same_dst_again) {
        const auto success = m_wf.m_store.lookup(dst_mbf, m_dst_walker);
        deterministic = success && m_dst_walker.m_deterministic.get(m_wf.iroot_part(ipart_dst));
        if (success) m_dst_walker.m_weight.copy_to(m_dst_weight);
        else m_dst_weight.assign(m_dst_weight.size(), dtype::null(m_dst_weight[0]));
    }
}

void Annihilator::annihilate_row(const uint_t &dst_ipart, const field::Mbf &dst_mbf, const wf_t &delta_weight,
                                 bool allow_initiation, Walker &dst_walker) {
    if (m_nadd == 0.0) {
        DEBUG_ASSERT_TRUE(allow_initiation,
                          "initiator rules are turned off, every initiating annihilation should be allowed");
    }
    DEBUG_ASSERT_FALSE(dst_mbf.is_zero(), "recv table is contiguous, shouldn't have any cleared rows");
    DEBUG_ASSERT_EQ(m_wf.m_dist.irank(dst_mbf), mpi::irank(),
                    "the received MBF has not been sent to the right MPI rank!")
    // zero magnitude weights should not have been communicated
    if (delta_weight == 0.0) return;

    m_wf.m_nspawned.m_local[dst_ipart] += std::abs(delta_weight);
    if (!dst_walker) {
        /*
         * the destination MBF row in m_wf.m_store is not currently occupied, so initiator rules must be applied
         */
        if (!allow_initiation) {
            //m_aborted_weight += std::abs(*delta_weight);
            return;
        }

        auto& new_walker = m_wf.create_row_(
                m_icycle, dst_mbf, m_prop.m_ham.get_energy(dst_mbf), m_refs.is_connected(dst_mbf));
        m_wf.set_weight(new_walker, dst_ipart, delta_weight);

    } else {
        wf_t weight_before = dst_walker.m_weight[dst_ipart];
        if (weight_before == 0.0 && !allow_initiation) {
            // the row exists, but that does not necessarily mean that each part has non-zero occupation
            //m_aborted_weight += std::abs(*delta_weight);
            return;
        }
        m_wf.m_nannihilated.m_local[dst_ipart] += annihilated_magnitude(weight_before, delta_weight);
        m_wf.change_weight(dst_walker, dst_ipart, delta_weight);
    }
}

void Annihilator::handle_dst_block(Spawn &block_begin, Spawn &next_block_begin,
                                   const wf_t &total_delta, Walker &dst_walker) {
    DEBUG_ASSERT_FALSE(in_same_dst_block(block_begin, next_block_begin),
                       "start of block and start of next block should not be in the same block");
    DEBUG_ASSERT_LT(block_begin.index(), next_block_begin.index(),
                    "block start should be strictly before current row in recv table");

    /*
     * only consider RDM contributions if there are any RDMs being accumulated and the destination exists
     */
    if (m_rdms && m_rdms.m_accum_epoch && dst_walker) {
        DEBUG_ASSERT_TRUE(block_begin.m_send_parents, "RDM sampling requires that parent MBFs are communicated");
        /*
         * store the original positions of the row objects in the recv table
         */
        const uint_t irow_begin = block_begin.index();
        const uint_t irow_next_begin = next_block_begin.index();

        auto &current = next_block_begin;
        wf_t src_weight = block_begin.m_src_weight;
        DEBUG_ONLY(src_weight);

        for (current.jump(irow_begin);; ++current) {
            if (!in_same_src_block(current, block_begin)) {
                handle_src_block(block_begin, dst_walker);
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
                DEBUG_ASSERT_EQ(src_weight, wf_t(current.m_src_weight),
                    "all spawns with the same dst_mbf, ipart_dst, and src_mbf should have exactly the same src_weight");
            }
        }
        /*
         * restore the original positions
         */
        block_begin.jump(irow_begin);
        next_block_begin.jump(irow_next_begin);
    }

    // spawning is never prevented if the destination walker exists and is populated in the relevant part
    bool allow_initiation = dst_walker && (dst_walker.m_weight[block_begin.m_ipart_dst] != 0.0);
    // always allow new walker to be created if multiple sources are simultaneously spawning here
    if (!allow_initiation) allow_initiation = block_begin.offset(next_block_begin) > 1;
    // else, only allow initiation if the lone parent was an initiator
    if (!allow_initiation) {
        // only one src_mbf for this dst_mbf. If the parent is an initiator,
        // contributions to unoccupied MBFs are allowed
        allow_initiation = block_begin.m_src_initiator;
    }
    annihilate_row(block_begin.m_ipart_dst, block_begin.m_dst_mbf, total_delta, allow_initiation, dst_walker);
    block_begin.jump(next_block_begin);
    DEBUG_ASSERT_EQ(next_block_begin.index(), block_begin.index(), "row not set to beginning of next block");
}

void Annihilator::handle_src_block(const Spawn &block_begin, const Walker &dst_row) {
    DEBUG_ASSERT_TRUE(m_rdms.m_accum_epoch, "shouldn't be sampling RDMs yet");
    DEBUG_ASSERT_EQ(block_begin.m_dst_mbf, dst_row.m_mbf, "wrong dst_row found");

    uint_t ipart_dst = block_begin.m_ipart_dst;

    /*
     * don't make contributions to RDM elements if they already take the equivalent contribution from deterministic
     * average connections to the Hartree-Fock determinant
     */
    if (m_hf) {
        if (dst_row.m_mbf == m_hf->mbf()) return;
        if (block_begin.m_src_mbf == m_hf->mbf()) return;
    }
    const auto iroot = m_wf.iroot_part(ipart_dst);
    /*
     * or if they already take contributions from deterministic subspace connections.
     * i.e. the src MBF is a deterministic subspace member for this root index, and so too is the dst MBF
     */
    if (block_begin.m_src_deterministic && dst_row.m_deterministic.get(iroot)) return;


    DEBUG_ASSERT_EQ(block_begin.m_dst_mbf, dst_row.m_mbf, "found row doesn't correspond to spawned dst");
    const auto ipart_replica = dst_row.ipart_replica(ipart_dst);
    wf_t contrib = m_dst_weight[ipart_replica];
    // recover pre-death value of replica population (on average)
    contrib /= 1.0 - m_prop.tau() * (dst_row.m_hdiag - m_prop.m_shift.m_values[ipart_replica]);
    contrib = arith::conj(contrib);
    contrib *= wf_t(block_begin.m_src_weight);
    m_rdms.make_contribs(block_begin.m_src_mbf, dst_row.m_mbf, contrib);
}

void Annihilator::loop_over_dst_mbfs() {
    if (!m_wf.recv().nrow_in_use()) return;
    /*
     * put the block_begin row to the beginning of the recv table
     */
    auto &block_begin = m_work_row1;
    block_begin.restart();

    bool dst_deterministic = false;
    wf_t total_delta = 0.0;
    m_dst_walker.select_null();

    /*
     * set m_dst_walker to the record corresponding to the child walker if it exists, and set the replica weights if the
     * current destination walker differs from the previous
     */
    lookup_dst(block_begin.m_dst_mbf, block_begin.m_ipart_dst, dst_deterministic);
    auto &current = m_work_row2;
    for (current.restart();; ++current) {
        if (!in_same_dst_block(current, block_begin)) {
            /*
             * either the current row has been iterated over the high-water mark of the recv table, or onto a row with a
             * different (dst_mbf, ipart_dst) pair. In either case, we have reached the end of a block of the major
             * sorting field, so we must handle the block just finished
             */
            handle_dst_block(block_begin, current, total_delta, m_dst_walker);
            DEBUG_ASSERT_EQ(block_begin.index(), current.index(),
                            "block_begin should have been pointed to the beginning of the next block");
            /*
             * the above call iterates block_begin through the contributions until the row index matches that of current
             * which would be the start of the next block if current is not past the end of recv rows. In any case, the
             * two rows should be the same. If the new start of block is past the end of recv rows, we're done
             */
            if (!block_begin) break;
            /*
             * new block of dst_mbfs, so reset total contrib
             */
            total_delta = 0.0;
            /*
             * and look it up in the wavefunction store
             * i.e. set m_dst_walker to the record corresponding to the child walker if it exists
             */
            lookup_dst(block_begin.m_dst_mbf, block_begin.m_ipart_dst, dst_deterministic);
        } else {
            if (current.m_send_parents) {
                DEBUG_ASSERT_NE(current.m_dst_mbf, current.m_src_mbf,
                                "should never have diagonal connections at annihilation");
            }
            DEBUG_ASSERT_EQ(current.m_dst_mbf, block_begin.m_dst_mbf, "dst MBFs should be the same");
            DEBUG_ASSERT_EQ(current.m_ipart_dst, block_begin.m_ipart_dst, "dst iparts should be the same");
        }
        if (!dst_deterministic || !current.m_src_deterministic) {
            // this is not a determ->determ connection, so include it
            total_delta += current.m_delta_weight;
        }
    }
}