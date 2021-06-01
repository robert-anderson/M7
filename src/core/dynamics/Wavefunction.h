//
// Created by Robert John Anderson on 2020-04-03.
//

#ifndef M7_WAVEFUNCTION_H
#define M7_WAVEFUNCTION_H


#include <src/core/io/Options.h>
#include <src/core/hamiltonian/Hamiltonian.h>
#include <src/core/parallel/Reduction.h>
#include "src/core/table/Communicator.h"
#include "src/core/dynamics/WalkerTable.h"
#include "src/core/dynamics/SpawnTable.h"
#include "src/core/parallel/RankAllocator.h"
#include "src/core/parallel/ReductionMember.h"
#include "src/core/field/Fields.h"
#include "src/core/sort/QuickSorter.h"

/**
 * A communicator whose m_store is the list of occupation number vectors currently
 * with a non-zero value in at least one of the elements of the NdNumberField m_weight.
 *
 * The wavefunction row WalkerTableRow contains multidimensional fields whose elements
 * are referred to as "parts".
 *
 * This m_store is often called the "main walker list", and the CommunicatingPair tables
 * are the called the "spawning lists"
 */
struct Wavefunction : Communicator<WalkerTableRow, SpawnTableRow> {

    const Options &m_opts;
    const size_t m_nsite;

    /**
     * all multidimensional array indices of the multidimensional fields of the m_store row
     */
    NdEnumeration<defs::ndim_wf> m_format;

    /**
     * collection of all reductions which are summed at the end of every cycle
     */
    ReductionSyndicate m_summables;

    /**
     * number of initiator ONVs in each part of the WF
     */
    NdReduction<size_t, defs::ndim_wf> m_ninitiator;
    /**
     * change over the last cycle in the number of initiator ONVs
     */
    NdReduction<int, defs::ndim_wf> m_delta_ninitiator;
    /**
     * number of ONVs with any associated weight in any part
     */
    Reduction<size_t> m_nocc_onv;
    /**
     * change in the number of occupied ONVs
     */
    Reduction<int> m_delta_nocc_onv;
    /**
     * L1 norm of each part of the WF
     */
    NdReduction<defs::wf_comp_t, defs::ndim_wf> m_nwalker;
    /**
     * change in the L1 norm
     */
    NdReduction<defs::wf_comp_t, defs::ndim_wf> m_delta_nwalker;
    /**
     * square of the L2 norm of each part of the WF
     */
    NdReduction<defs::wf_comp_t, defs::ndim_wf> m_l2_norm_square;
    /**
     * change in the L2 norm
     */
    NdReduction<defs::wf_comp_t, defs::ndim_wf> m_delta_l2_norm_square;
    /**
     * number of walkers received in spawning process
     */
    NdReduction<defs::wf_comp_t, defs::ndim_wf> m_nspawned;
    /**
     * number of walkers annihilated in the loop_over_spawned method for each part
     */
    NdReduction<defs::wf_comp_t, defs::ndim_wf> m_nannihilated;

    Wavefunction(const Options &opts, size_t nsite);

    std::vector<std::string> h5_field_names();

    void h5_write(hdf5::GroupWriter &parent, std::string name = "wavefunction");

    void h5_read(hdf5::GroupReader &parent, const Hamiltonian<> &ham, const fields::Onv<> &ref,
                 std::string name = "wavefunction");


    void begin_cycle();

    void end_cycle();

//    void update(size_t icycle, double work_time) {
//        m_ra.update(icycle, work_time, m_walkers, m_walkers.m_key_field);
//    }

    const size_t& nroot() const {
        return m_format.extent(0);
    }

    const size_t& nreplica() const {
        return m_format.extent(1);
    }

    size_t ipart_replica(const size_t& ipart) const {
        return nreplica()==1 ? ipart : (ipart/2)*2+!(ipart&1ul);
    }

    defs::wf_comp_t square_norm(const size_t& ipart) const;

    defs::wf_comp_t l1_norm(const size_t& ipart) const;

    /**
     * allow the current ONV in m_store.m_row to change the weight on ONVs to which it generates
     * spawned contributions, and update initiator statistics accordingly
     * @param ipart
     *  flat index of the initiator flag to be set in the selected row
     */
    void grant_initiator_status(const size_t &ipart);

    /**
     * prohibit the current ONV in m_store.m_row to change the weight on ONVs to which it generates
     * spawned contributions, and update initiator statistics accordingly
     * @param ipart
     *  flat index of the initiator flag to be cleared in the selected row
     */
    void revoke_initiator_status(const size_t &ipart);

    /**
     * all changes in the m_weight member of any row associated with m_store should occur through
     * this function so that changes can be properly recorded
     * @param ipart
     *  part index of the WF to update
     * @param new_weight
     *  value to which this part weight is to be set
     */
    void set_weight(const size_t &ipart, const defs::wf_t &new_weight);

    void set_weight(const defs::wf_t &new_weight) {
        for (size_t ipart=0ul; ipart<m_format.nelement(); ++ipart) set_weight(ipart, new_weight);
    }

    void set_weight(const fields::Numbers<defs::wf_t, defs::ndim_wf> &new_weight){
        for (size_t i=0ul; i < m_format.nelement(); ++i) set_weight(i, new_weight[i]);
    }

    /**
     * convenience method to set_weight based on a difference relative to the current weight of
     * the part
     * @param ipart
     *  part index
     * @param delta
     *  change in the weight
     */
    void change_weight(const size_t &ipart, const defs::wf_t &delta);

    /**
     * convenience method to set_weight based on a scalar factor relative to current weight
     * @param ipart
     *  part index
     * @param factor
     *  fractional change in the weight
     */
    void scale_weight(const size_t &ipart, const double &factor);

    /**
     * convenience method to set the weight of a part of the WF to zero on the currently
     * selected row
     * @param ipart
     *  part index
     */
    void zero_weight(const size_t &ipart);

    void remove_row();

private:
    /**
     * Only called on the rank assigned to the ONV by the RankAllocator
     * @param icycle
     *  MC cycle index on which ONV is being added
     * @param onv
     *  ONV of row to be added
     * @param hdiag
     *  diagonal matrix element is cached here
     * @return
     *  index of created row
     * @param refconn
     *  element true if reference ONV of corresponding WF part is connected
     *  i.e. the connection to the reference has a non-zero H matrix element
     * @return
     */
    size_t create_row_(const size_t &icycle, const fields::Onv<> &onv,
                       const defs::ham_comp_t &hdiag, const std::vector<bool>& refconns) {
        ASSERT(refconns.size()==npart());
        ASSERT(mpi::i_am(m_ra.get_rank(onv)));
        if (m_store.is_full()) m_store.expand(1);
        auto irow = m_store.insert(onv);
        m_delta_nocc_onv.m_local++;
        m_store.m_row.jump(irow);
        ASSERT(m_store.m_row.m_onv == onv)
        m_store.m_row.m_hdiag = hdiag;
        for (size_t ipart=0ul; ipart<npart(); ++ipart)
            m_store.m_row.m_reference_connection.put(ipart, refconns[ipart]);
        if (defs::enable_mevs) {
            m_store.m_row.m_icycle_occ = icycle;
            m_store.m_row.m_average_weight = 0;
        }
        return irow;
    }


    size_t create_row_(const size_t &icycle, const fields::Onv<> &onv,
                       const defs::ham_comp_t &hdiag, bool refconn) {
        return create_row_(icycle, onv, hdiag, std::vector<bool>(npart(), refconn));
    }

public:
    /**
     * Called on all ranks, dispatching create_row_ on the assigned rank only
     */
    TableBase::Loc create_row(const size_t& icycle, const fields::Onv<> &onv,
                              const defs::ham_comp_t &hdiag, const std::vector<bool>& refconns) {
        size_t irank = m_ra.get_rank(onv);
        size_t irow;
        if (mpi::i_am(irank)) irow = create_row_(icycle, onv, hdiag, refconns);
        mpi::bcast(irow, irank);
        return {irank, irow};
    }


    TableBase::Loc create_row(const size_t& icycle, const fields::Onv<> &onv,
                              const defs::ham_comp_t &hdiag, bool refconn) {
        return create_row(icycle, onv, hdiag, std::vector<bool>(npart(), refconn));
    }

    size_t add_spawn(const fields::Onv<> &dst_onv, const defs::wf_t &delta,
                     bool initiator, bool deterministic, size_t dst_ipart);

    size_t add_spawn(const fields::Onv<> &dst_onv, const defs::wf_t &delta,
                     bool initiator, bool deterministic, size_t dst_ipart,
                     const fields::Onv<> &src_onv, const defs::wf_t &src_weight);

    const size_t& npart() const {
        return m_format.nelement();
    }

};

#endif //M7_WAVEFUNCTION_H
