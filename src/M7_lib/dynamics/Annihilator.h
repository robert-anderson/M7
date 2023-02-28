//
// Created by Robert J. Anderson on 14/07/2021.
//

#ifndef M7_ANNIHILATOR_H
#define M7_ANNIHILATOR_H

#include <M7_lib/wavefunction/Wavefunction.h>
#include <M7_lib/wavefunction/Reference.h>
#include <M7_lib/bilinear/Rdms.h>

/**
 * updates the state of the stored wavefunction subject to the new walkers in the receive buffer,
 * also sorts communicated walkers according to their destination many-body basis function and makes contributions to
 * MAEs sue to the sampled connections
 */
struct Annihilator {
    /**
     * wavefunction whose recv and store tables are merged herein
     */
    wf::Fci &m_wf;
    /**
     * propagator object providing access to Hamiltonian and deterministic subspace information
     */
    const Propagator& m_prop;
    /**
     * references are needed to discern whether newly-created walkers contribute to the reference energies
     */
    const wf::Refs& m_refs;
    /**
     * Hartree-Fock MBF, fixed throughout the calculation, for Brillouin theorem RDM accumulation, may be null
     */
    const shared_rows::Walker* m_hf;
    /**
     * all RDM objects, including contracted tensors
     */
    Rdms& m_rdms;
    /**
     * initiator threshold
     */
    const wf_comp_t m_nadd;
    /**
     * reference to the current cycle counter of the Solver class
     */
    const uint_t& m_icycle;

private:
    /**
     * traversal of the wavefunction recv table requires two working rows, for the beginning of adjacent blocks
     * a pair of working rows is also needed for sorting the received spawns
     */
    Spawn m_work_row1, m_work_row2;
    /**
     * Row pointing to the destination walker record if it is found to exist
     */
    Walker m_dst_walker;
    /**
     * weights of the dst walker are cached so their pre-annihilation values are available to MAE accumulation
     */
    v_t<wf_t> m_dst_weight;
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
        DEBUG_ASSERT_TRUE(row1 || row2, "rows should not both be out of dereferencable range");
        if (bool(row1) != bool(row2)) return false;
        return (row1.m_dst_mbf == row2.m_dst_mbf) && (row1.m_ipart_dst == row2.m_ipart_dst);
    }

    static bool in_same_src_block(const Spawn& row1, const Spawn& row2){
        if (!in_same_dst_block(row1, row2)) return false;
        return row1.m_src_mbf == row2.m_src_mbf;
    }

public:

    Annihilator(wf::Fci &wf, const Propagator& prop, const wf::Refs& refs,
                const shared_rows::Walker* hf, Rdms& rdms, const uint_t& icycle, wf_comp_t nadd);

    /**
     * using the comparator implementation, sort the m_recv table of the referenced wavefunction
     */
    void sort_recv();

    /**
     * lookup the destination walker row in the wavefunction store table and point m_dst_walker at the result.
     * also, cache the weights on all parts of the dst walker so sequential annihilation of replicas can take place
     * without biasing RDM contributions
     * @param dst_mbf
     *  destination many-body basis function
     * @param ipart_dst
     *  index among walker populations for which the spawned contribution is destined
     * @param deterministic
     *  return with value true if the dst walker exists and is in the deterministic subspace
     */
    void lookup_dst(const Mbf& dst_mbf, uint_t ipart_dst, bool& deterministic);

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

    void handle_src_block(const Spawn &block_begin, const Walker& dst_row);

    /**
     * loop through all received spawned rows and group them into blocks using a pair of Row objects. Once a block has
     * been identified, handle its MAE and WF contributions via handle_dst_block method
     */
    void loop_over_dst_mbfs();

    /**
     * TODO: complex sign problem or "phase problem"
     */
    template<typename T>
    static T annihilated_magnitude(const T& weight, const T& delta){
        if (weight > 0.0 && delta < 0.0) return std::min(weight, -delta);
        if (weight < 0.0 && delta > 0.0) return std::min(-weight, delta);
        return 0.0;
    }

    template<typename T>
    static T annihilated_magnitude(const std::complex<T>& weight, const std::complex<T>& delta){
        return annihilated_magnitude(arith::real_ref(weight), arith::real_ref(delta)) +
               annihilated_magnitude(arith::imag_ref(weight), arith::imag_ref(delta));
    }
};


#endif //M7_ANNIHILATOR_H
