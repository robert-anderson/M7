//
// Created by rja on 14/07/2021.
//

#include "Annihilator.h"

comparators::index_cmp_fn_t Annihilator::make_sort_cmp_fn() {
    if (m_send_parents) {
        return [&](const size_t &irow1, const size_t &irow2) {
            m_work_row1.jump(irow1);
            m_work_row2.jump(irow2);
            // sort criteria from major to minor: dst ONV, dst_ipart, src ONV,
            if (m_work_row1.m_dst_onv == m_work_row2.m_dst_onv) {
                if (m_work_row1.m_dst_ipart == m_work_row2.m_dst_ipart) {
                    return m_work_row1.m_src_onv < m_work_row2.m_src_onv;
                }
                return m_work_row1.m_dst_ipart < m_work_row2.m_dst_ipart;
            }
            return m_work_row1.m_dst_onv < m_work_row2.m_dst_onv;
        };
    } else {
        return [&](const size_t &irow1, const size_t &irow2) {
            m_work_row1.jump(irow1);
            m_work_row2.jump(irow2);
            // sort criteria from major to minor: dst ONV, dst_ipart, src ONV,
            if (m_work_row1.m_dst_onv == m_work_row2.m_dst_onv) {
                return m_work_row1.m_dst_ipart < m_work_row2.m_dst_ipart;
            }
            return m_work_row1.m_dst_onv < m_work_row2.m_dst_onv;
        };
    }
}

Annihilator::Annihilator(Wavefunction &wf, const Hamiltonian<> &ham, const References &refs,
                         const size_t &icycle, defs::wf_t nadd) :
        m_wf(wf), m_ham(ham), m_refs(refs), m_nadd(nadd), m_icycle(icycle),
        m_work_row1(wf.m_comm.recv().m_row), m_work_row2(wf.m_comm.recv().m_row),
        m_send_parents(m_work_row1.m_send_parents), m_sort_cmp_fn(make_sort_cmp_fn()) {}

void Annihilator::sort_recv() {
    QuickSorter qs(m_sort_cmp_fn);
    qs.reorder_sort(m_wf.recv());
}

void
Annihilator::annihilate_row(const size_t &dst_ipart, const fields::Onv<> &dst_onv, const defs::wf_t &delta_weight,
                            bool allow_initiation, const size_t &irow_store) {
    if (m_nadd == 0.0) {
        DEBUG_ASSERT_TRUE(allow_initiation,
                          "initiator rules are turned off, every initiating annihilation should be allowed");
    }
    DEBUG_ASSERT_FALSE(dst_onv.is_zero(), "recv table is contiguous, shouldn't have any cleared rows");
    DEBUG_ASSERT_EQ(m_wf.m_ra.get_rank(dst_onv), mpi::irank(),
                    "the received ONV has not been sent to the right MPI rank!")
    // zero magnitude weights should not have been communicated
    if (consts::float_is_zero(delta_weight)) return;

    m_wf.m_nspawned.m_local[dst_ipart] += std::abs(delta_weight);
    if (irow_store == ~0ul) {
        /*
         * the destination ONV is not currently occupied, so initiator rules
         * must be applied
         */
        if (!allow_initiation) {
            //m_aborted_weight += std::abs(*delta_weight);
            return;
        }

        m_wf.create_row_(m_icycle, dst_onv, m_ham.get_energy(dst_onv), m_refs.is_connected(dst_onv));
        m_wf.set_weight(dst_ipart, delta_weight);

    } else {
        m_wf.m_store.m_row.jump(irow_store);
        defs::wf_t weight_before = m_wf.m_store.m_row.m_weight[dst_ipart];
        auto weight_after = weight_before + delta_weight;
        if (!consts::float_is_zero(weight_before) && !consts::float_is_zero(weight_after)
            && ((std::abs(weight_before) > 0) != (std::abs(weight_after) > 0)))
            m_wf.m_nannihilated.m_local[dst_ipart] += std::abs(std::abs(weight_before) - std::abs(weight_after));
        m_wf.change_weight(dst_ipart, delta_weight);
    }
}

void
Annihilator::handle_dst_onv_block(SpawnTableRow &block_start, SpawnTableRow &current, const defs::wf_t &total_delta,
                                  const size_t &irow_store) {
    DEBUG_ASSERT_LT(block_start.index(), current.index(),
                    "block start should be strictly before current row in recv table");
    // row_block_start is now at last row in last block
    bool allow_initiation = (current.index() - block_start.index()) > 1;
    if (!allow_initiation) {
        // only one src_onv for this dst_onv. If the parent is an initiator,
        // contributions to unoccupied ONVs are allowed
        allow_initiation = block_start.m_src_initiator;
    }
    annihilate_row(0, block_start.m_dst_onv, total_delta, allow_initiation, irow_store);
    block_start.jump(current);
}

void Annihilator::loop_over_dst_onvs() {
    if (!m_wf.recv().m_hwm) return;

    auto &block_start = m_work_row1;
    block_start.restart();
    defs::wf_t total_delta = 0.0;

    auto &store_row = m_wf.m_store.m_row;
    bool dst_deterministic = false;
    size_t irow_dst = ~0ul;
    auto find_dst = [&]() {
        irow_dst = *m_wf.m_store[block_start.m_dst_onv];
        if (irow_dst!=~0ul) {
            store_row.jump(irow_dst);
            dst_deterministic = store_row.m_deterministic.get(0);
        }
        else dst_deterministic = false;
    };

    find_dst();

    auto &current = m_work_row2;
    for (current.restart(); ;current.step()) {
        if (!current.in_range() || (current.m_dst_onv != block_start.m_dst_onv)) {
            // handle the block just finished
            handle_dst_onv_block(block_start, current, total_delta, irow_dst);
            DEBUG_ASSERT_EQ(block_start.index(), current.index(),
                            "block_start should have been pointed to the beginning of the next block");
            if (!block_start.in_range()) return;
            total_delta = 0.0;
            find_dst();
        } else {
            DEBUG_ASSERT_EQ(current.m_dst_onv, block_start.m_dst_onv, "dst ONVs should be the same");
        }
        if (!dst_deterministic || !current.m_src_deterministic) {
            // this is not a determ->determ connection, so include it
            total_delta += current.m_delta_weight;
        }
    }
}
