//
// Created by Robert J. Anderson on 14/07/2021.
//

#ifndef M7_ANNIHILATOR_H
#define M7_ANNIHILATOR_H


#include <M7_lib/wavefunction/Wavefunction.h>
#include <M7_lib/wavefunction/Reference.h>
#include <M7_lib/bilinear/Rdm.h>

/**
 * we need to be able to lookup the dst in the main table.
 */
struct DstLookup {
    /**
     * wavefunction object
     */
    Wavefunction& m_wf;
    /**
     * true if destination Walker is in the deterministic subspace
     */
    bool m_deterministic = false;
    /**
     * row object for pointing to the destination Walker record (or null if not found)
     */
    Walker m_row;
    DstLookup(Wavefunction& wf);

    /**
     * @param m_block_start_row
     *  row pointing to the first record in the current block in the recv table
     * @return
     */
    Lookup lookup(const Spawn &recv_row);
};


/**
 * updates the state of the stored wavefunction subject to the new walkers in the receive buffer,
 * also sorts communicated walkers according to their destination many-body basis function and makes contributions to
 * MAEs sue to the sampled connections
 */
struct Annihilator {
    Wavefunction &m_wf;
    const Propagator& m_prop;
    const References& m_refs;
    Rdms& m_rdms;
    const wf_comp_t m_nadd;
    const uint_t& m_icycle;

private:
    Spawn m_work_row1;
    Spawn m_work_row2;
    Walker m_dst_walker;
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

    /**
     * @param row1
     *  spawning row
     * @param row2
     *  another spawning row to compare with
     * @return
     *  true if row1 and row2 are in distinct (dst_mbf, ipart_dst) blocks
     */
    static bool in_same_dst_block(const Spawn& row1, const Spawn& row2){
        DEBUG_ASSERT_TRUE(row1.in_range() || row2.in_range(), "rows should not both be out of range");
        if (row1.in_range() != row2.in_range()) return false;
        return (row1.m_dst_mbf == row2.m_dst_mbf) && (row1.m_ipart_dst == row2.m_ipart_dst);
    }

    static bool in_same_src_block(const Spawn& row1, const Spawn& row2){
        if (!in_same_dst_block(row1, row2)) return false;
        return row1.m_src_mbf == row2.m_src_mbf;
    }

public:

    Annihilator(Wavefunction &wf, const Propagator& prop, const References& refs, Rdms& rdms,
                const uint_t& icycle, wf_comp_t nadd);

    /**
     * using the comparator implementation, sort the m_recv table of the referenced wavefunction
     */
    void sort_recv();

    /**
     * lookup the destination walker row in the wavefunction store table and store the result in m_dst_walker
     * @param recv_row
     *  row in the spawning recv table
     * @param deterministic
     *  return with value true if the dst walker exists and is in the deterministic subspace
     * @return
     *  lookup result
     */
    Lookup lookup_dst(const Spawn& recv_row, bool& deterministic);

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
     */
    void annihilate_row(const uint_t &dst_ipart, const field::Mbf &dst_mbf, const wf_t &delta_weight,
                        bool allow_initiation, Walker& dst_walker);
    /**
     * given that all rows between block_start (inclusively) and current (exclusively) correspond the the same (dst_mbf,
     * ipart_dst) pair, make all contributions to that pair in m_wf.m_store only after making all contributions to the
     * MAEs. block_start should be pointing to the same row as current on termination of this function
     * @param block_begin
     *  first row of the (dst_mbf, ipart_dst) block being handled
     * @param next_block_begin
     *  first row of the next (dst_mbf, ipart_dst) block or past the end of the recv table
     * @param total_delta
     *  total change in the ipart_dst indexed weight element of dst_mbf due to the sum of all spawns
     * @param dst_walker
     *  row in m_wf.m_store which stores the dst_mbf if found
     */
    void handle_dst_block(Spawn &block_begin, Spawn &next_block_begin,
                          const wf_t &total_delta, Walker& dst_walker);

    void handle_src_block(Spawn &block_begin, Walker& dst_row);

    /**
     * loop through all received spawned rows and group them into blocks using a pair of Row objects. Once a block has
     * been identified, handle its MAE and WF contributions via handle_dst_block method
     */
    void loop_over_dst_mbfs();

    /**
     * TODO: complex sign problem or "phase problem"
     *
     *  magnitude of the sign difference between two real-valued walker weights
     */
    template<typename T>
    static T annihilated_magnitude(const T& weight, const T& delta){
        if (weight > 0.0 && delta < 0.0) return std::min(weight, -delta);
        if (weight < 0.0 && delta > 0.0) return std::min(-weight, delta);
        return 0.0;
    }
};


#endif //M7_ANNIHILATOR_H
