//
// Created by Robert John Anderson on 2020-04-03.
//

#ifndef M7_WAVEFUNCTION_H
#define M7_WAVEFUNCTION_H

#include <src/core/io/Options.h>
#include <src/core/hamiltonian/Hamiltonian.h>
#include <src/core/parallel/Reduction.h>
#include <src/core/io/Archivable.h>
#include "src/core/table/Communicator.h"
#include "src/core/dynamics/WalkerTable.h"
#include "src/core/dynamics/SpawnTable.h"
#include "src/core/parallel/RankAllocator.h"
#include "src/core/parallel/ReductionMember.h"
#include "src/core/field/Fields.h"
#include "src/core/sort/QuickSorter.h"
#include "src/core/sort/GlobalExtremalRows.h"

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
struct Wavefunction : Communicator<WalkerTableRow, SpawnTableRow>, Archivable {

    typedef GlobalExtremalRows<WalkerTableRow, defs::wf_t, defs::ndim_wf> weights_gxr_t;

    const fciqmc_config::Document &m_opts;
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

    Wavefunction(const fciqmc_config::Document &opts, size_t nsite);

    ~Wavefunction();

    static bool need_send_parents(const fciqmc_config::Document &opts) {
        return opts.m_av_ests.m_fermion_rdm.m_rank > 0;
    }

    static bool need_av_weights(const fciqmc_config::Document &opts) {
        if (opts.m_av_ests.m_fermion_rdm.m_rank > 0) return true;
        return opts.m_av_ests.m_ref_excits.m_max_exlvl > 0;
    }

    bool storing_av_weights() const {
        return m_store.m_row.m_average_weight.belongs_to_row();
    }

    std::vector<std::string> h5_field_names();

    void h5_write(hdf5::GroupWriter &parent, std::string name = "wavefunction");

    void h5_read(hdf5::GroupReader &parent, const Hamiltonian<> &ham, const fields::Onv<> &ref,
                 std::string name = "wavefunction");

    void begin_cycle();

    void end_cycle();

    const size_t& nroot() const {
        return m_format.m_shape[0];
    }

    const size_t& nreplica() const {
        return m_format.m_shape[1];
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
        for (size_t ipart=0ul; ipart<m_format.m_nelement; ++ipart) set_weight(ipart, new_weight);
    }

    void set_weight(const fields::Numbers<defs::wf_t, defs::ndim_wf> &new_weight){
        for (size_t i=0ul; i < m_format.m_nelement; ++i) set_weight(i, new_weight[i]);
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

    void sort_recv();

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
            m_store.m_row.m_ref_conn.put(ipart, refconns[ipart]);
        if (storing_av_weights()) {
            m_store.m_row.m_icycle_occ = icycle;
            m_store.m_row.m_average_weight = 0;
        }
        return irow;
    }


    size_t create_row_(const size_t &icycle, const fields::Onv<> &onv,
                       const defs::ham_comp_t &hdiag, bool refconn) {
        return create_row_(icycle, onv, hdiag, std::vector<bool>(npart(), refconn));
    }

    /**
     * Called on all ranks, dispatching create_row_ on the assigned rank only
     */
    TableBase::Loc create_row(const size_t& icycle, const fields::Onv<> &onv,
                              const defs::ham_comp_t &hdiag, const std::vector<bool>& refconns) {
        size_t irank = m_ra.get_rank(onv);
        size_t irow;
        if (mpi::i_am(irank)) {
            irow = create_row_(icycle, onv, hdiag, refconns);
        }
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
        return m_format.m_nelement;
    }

private:

    void orthogonalize(NdReduction<defs::wf_t, 3>& overlaps,
                       const size_t& iroot, const size_t& jroot, const size_t& ireplica) {
        ASSERT(iroot<=jroot);
        auto& row = m_store.m_row;
        const auto ipart_src = m_format.flatten({iroot, ireplica});
        const auto ipart_dst = m_format.flatten({jroot, ireplica});
        overlaps.m_local[{iroot, jroot, ireplica}] +=
                consts::conj(row.m_weight[ipart_src])*row.m_weight[ipart_dst];
        if (jroot+1<nroot()) {
            // there is another part to project onto
            const auto ipart_next = m_format.flatten({jroot+1, ireplica});
            overlaps.m_local[{iroot, jroot + 1, ireplica}] +=
                    consts::conj(row.m_weight[ipart_src]) * row.m_weight[ipart_next];
        }
        if (iroot<jroot){
            const auto& overlap = overlaps.m_reduced[{iroot, jroot, ireplica}];
            const auto& norm = overlaps.m_reduced[{iroot, iroot, ireplica}];
            ASSERT(std::abs(norm)>1e-12);
            const auto gs_coeff = overlap / norm;
            change_weight(ipart_dst, -gs_coeff*row.m_weight[ipart_src]);
        }
    }
public:

    /**
     * using the Gram-Schmidt method, orthogonalize each (replicated) population with respect to all lower-lying states
     * by projecting their orthogonal complements.
     * the GS process transforms a non-orthogonal set {v} into an orthogonal set {u} by the algorithm:
     *  u0 = v0
     *  u1 = v1 - p(u0, v1)
     *  u2 = v2 - p(u0, v2) - p(u1, v2)
     *  .
     *  .
     *  .
     *  uk = vk - sum_{i=0}^{k-1} p(ui, vk)
     *
     *  where p(u, k) is the projection (<u, k>/<u, u> u)
     *
     * when calling this function, it is assumed that the value of m_l2_norm_square corresponds to the current walker
     * list.
     *
     * for i in [0, nroot)
     *  project psi_0 out of all higher states
     *
     */
    void orthogonalize() {
        // bra root, ket root, replica
        NdReduction<defs::wf_t, 3> overlaps({{nroot(), nroot(), nreplica()}});
        auto& row = m_store.m_row;
        for (size_t iroot = 0ul; iroot < nroot(); ++iroot) {
            for (size_t jroot = iroot; jroot < nroot(); ++jroot) {
                for (size_t ireplica = 0ul; ireplica < nreplica(); ++ireplica) {
                    for (row.restart(); row.in_range(); row.step()) {
                        if (!row.m_onv.is_zero()) orthogonalize(overlaps, iroot, jroot, ireplica);
                    }
                    overlaps.all_sum();
                }
            }
        }
    }


protected:
    void load_fn(hdf5::GroupReader &parent) override;

    void save_fn(hdf5::GroupWriter &parent) override;

};

#endif //M7_WAVEFUNCTION_H
