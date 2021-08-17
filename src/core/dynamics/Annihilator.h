//
// Created by rja on 14/07/2021.
//

#ifndef M7_ANNIHILATOR_H
#define M7_ANNIHILATOR_H


#include <src/core/wavefunction/Wavefunction.h>
#include <src/core/wavefunction/Reference.h>

/**
 * we need to be able to lookup the dst in the main table.
 */
struct DstFinder {
    Wavefunction& m_wf;
    WalkerTableRow &m_store_row;
    SpawnTableRow &m_block_start_row;
    bool m_dst_deterministic = false;
    size_t m_irow_dst = ~0ul;
    DstFinder(Wavefunction& wf, SpawnTableRow& block_start_row):
    m_store_row(wf.m_store.m_row), m_block_start_row(block_start_row){}

    void operator()() {
        irow_dst = *m_wf.m_store[block_start.m_dst_mbf];
        if (irow_dst != ~0ul) {
            store_row.jump(irow_dst);
            dst_deterministic = store_row.m_deterministic.get(0);
        } else dst_deterministic = false;
    };

};


/**
 * updates the state of the stored wavefunction subject to the new walkers in the receive buffer,
 * also sorts communicated walkers according to their destination many-body basis function and makes contributions to
 * MAEs sue to the sampled connections
 */
struct Annihilator {
    Wavefunction &m_wf;
    const Hamiltonian& m_ham;
    const References& m_refs;
    const defs::wf_comp_t m_nadd;
    const size_t& m_icycle;

private:
    SpawnTableRow m_work_row1;
    SpawnTableRow m_work_row2;
    const bool m_send_parents;
    /**
     * function that returns the relative order of two rows in the receive buffer
     */
    const comparators::index_cmp_fn_t m_sort_cmp_fn;
    /**
     * the comparator function has two different implementations: one which sorts only according to the destination MBF
     * and the WF part which takes the contribution, and another implementation which additionally sorts those blocks
     * by the source (or "parent") MBF
     * @return
     *  index comparator function
     */
    comparators::index_cmp_fn_t make_sort_cmp_fn();

public:

    Annihilator(Wavefunction &wf, const Hamiltonian& ham, const References& refs, const size_t& icycle, defs::wf_comp_t nadd);

    /**
     * using the comparator implementation, sort the m_recv table of the referenced wavefunction
     */
    void sort_recv();

    /**
     * performs coefficient lookup and (if initiator criteria are satisfied) update due to the sum of all spawned
     * contributions to the same dst_mbf, dst_ipart pair
     * @param dst_ipart
     *  part index of the wavefunction due to receive the total contribution
     * @param dst_mbf
     *  many-body basis function object due to receive the total contribution
     * @param delta_weight
     *  total contribution for all occurrences of the dst_mbf, dst_ipart pair in the recv table
     * @param allow_initiation
     *  if the lookup of the dst_mbf does not find an existing row in the m_wf.m_store, or row does exist but the weight
     *  at index dst_ipart is zero, then the contribution should only be made if this arg is true, otherwise the spawn
     *  is said to be aborted. Allow initiation is true if any of the contributing src_mbfs are initiators, or if there
     *  are multiple src_mbfs spawning to the same zero-weighted part in the same cycle
     * @param irow_store
     *  the result of the lookup on the main m_wf.m_store table
     */
    void annihilate_row(const size_t &dst_ipart, const field::Mbf &dst_mbf, const defs::wf_t &delta_weight,
                        bool allow_initiation, const size_t &irow_store);

    void handle_dst_mbf_block(SpawnTableRow &block_start, SpawnTableRow &current,
                              const defs::wf_t &total_delta, const size_t& irow_store);

    void loop_over_dst_mbfs();

};


#endif //M7_ANNIHILATOR_H
