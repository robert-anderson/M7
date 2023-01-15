//
// Created by Robert John Anderson on 2020-04-03.
//

#ifndef M7_WAVEFUNCTION_H
#define M7_WAVEFUNCTION_H

#include <M7_lib/io/Options.h>
#include <M7_lib/hamiltonian/Hamiltonian.h>
#include <M7_lib/parallel/Reduction.h>
#include <M7_lib/io/Archivable.h>
#include <M7_lib/communication/Communicator.h>
#include <M7_lib/parallel/RankAllocator.h>
#include <M7_lib/field/Fields.h>
#include <M7_lib/sort/QuickSort.h>
#include <M7_lib/sort/GlobalExtremalRows.h>

#include "WalkerTable.h"
#include "SpawnTable.h"
#include "FciInitializer.h"

/**
 * A communicator whose m_store is the list of occupation number vectors currently
 * with a non-zero value in at least one of the elements of the NdNumberField m_weight.
 *
 * The wavefunction row Walker contains multidimensional fields whose elements
 * are referred to as "parts".
 *
 * This m_store is often called the "main walker list", and the SendRecv tables
 * are the called the "spawning lists"
 */
struct Wavefunction : communicator::BasicSend<Walker, Spawn> {
    typedef GlobalExtremalRows<wf_t, c_ndim_wf> weights_gxr_t;

    const conf::Document &m_opts;
    const sys::Sector m_sector;

    /**
     * all multidimensional array indices of the multidimensional fields of the m_store row
     */
    NdEnumeration<c_ndim_wf> m_format;

    /**
     * collection of all reductions which are summed at the end of every cycle
     */
    ReductionSyndicate m_summables;

    /**
     * number of initiator MBFs in each part of the WF
     */
    NdReduction<uint_t, c_ndim_wf> m_ninitiator;
    /**
     * number of MBFs with any associated weight in any part
     */
    Reduction<uint_t> m_nocc_mbf;
    /**
     * change in the number of occupied MBFs
     */
    Reduction<int> m_delta_nocc_mbf;
    /**
     * L1 norm of each part of the WF
     */
    NdReduction<wf_comp_t, c_ndim_wf> m_nwalker;
    /**
     * change in the L1 norm
     */
    NdReduction<wf_comp_t, c_ndim_wf> m_delta_nwalker;
    /**
     * square of the L2 norm of each part of the WF
     */
    NdReduction<wf_comp_t, c_ndim_wf> m_l2_norm_square;
    /**
     * change in the L2 norm
     */
    NdReduction<wf_comp_t, c_ndim_wf> m_delta_l2_norm_square;
    /**
     * number of walkers received in spawning process
     */
    NdReduction<wf_comp_t, c_ndim_wf> m_nspawned;
    /**
     * number of walkers annihilated in the loop_over_spawned method for each part
     */
    NdReduction<wf_comp_t, c_ndim_wf> m_nannihilated;

    Wavefunction(const conf::Document &opts, const sys::Sector& sector);

    ~Wavefunction();

    void log_top_weighted(uint_t ipart, uint_t nrow=20);

    static bool need_send_parents(const conf::Document &opts) {
        return opts.m_av_ests.m_rdm.m_ranks.m_value.size();
    }

    static bool need_av_weights(const conf::Document &opts) {
        if (need_send_parents(opts)) return true;
        return opts.m_av_ests.m_hf_excits.m_max_exlvl > 0;
    }

    bool storing_av_weights() const {
        return m_store.m_row.m_average_weight.belongs_to_row();
    }

    void h5_write(const hdf5::NodeWriter &parent, str_t name = "wavefunction");

    void h5_read(const hdf5::NodeReader &parent, const Hamiltonian &ham, const field::Mbf &ref,
                 str_t name = "wavefunction");

    void begin_cycle();

    void end_cycle();

    uint_t nroot() const {
        return m_store.m_row.nroot();
    }

    uint_t nreplica() const {
        return m_store.m_row.nreplica();
    }

    uint_t ipart_replica(uint_t ipart) const {
        return m_store.m_row.ipart_replica(ipart);
    }

    uint_t iroot_part(uint_t ipart) const {
        return ipart/2;
    }

    wf_comp_t square_norm(uint_t ipart) const;

    wf_comp_t l1_norm(uint_t ipart) const;

    /**
     * all changes in the m_weight member of any row associated with m_store should occur through
     * this function so that changes can be properly recorded
     * @param ipart
     *  part index of the WF to update
     * @param new_weight
     *  value to which this part weight is to be set
     */
    void set_weight(Walker& walker, uint_t ipart, wf_t new_weight);

    void set_weight(Walker& walker, wf_t new_weight) {
        for (uint_t ipart=0ul; ipart<m_format.m_nelement; ++ipart) set_weight(walker, ipart, new_weight);
    }

    void set_weight(Walker& walker, const field::Numbers<wf_t, c_ndim_wf> &new_weight){
        for (uint_t i=0ul; i < m_format.m_nelement; ++i) set_weight(walker, i, new_weight[i]);
    }

    /**
     * convenience method to set_weight based on a difference relative to the current weight of
     * the part
     * @param ipart
     *  part index
     * @param delta
     *  change in the weight
     */
    void change_weight(Walker& walker, uint_t ipart, wf_t delta);

    /**
     * convenience method to set_weight based on a scalar factor relative to current weight
     * @param ipart
     *  part index
     * @param factor
     *  fractional change in the weight
     */
    void scale_weight(Walker& walker, uint_t ipart, double factor);

    /**
     * convenience method to set the weight of a part of the WF to zero on the currently
     * selected row
     * @param ipart
     *  part index
     */
    void zero_weight(Walker& walker, uint_t ipart);

    void remove_row(Walker& walker);

    /**
     * Only called on the rank assigned to the MBF by the RankAllocator
     * @param icycle
     *  MC cycle index on which MBF is being added
     * @param mbf
     *  MBF of row to be added
     * @param hdiag
     *  diagonal matrix element is cached here
     * @return
     *  index of created row
     * @param refconn
     *  element true if reference MBF of corresponding WF part is connected
     *  i.e. the connection to the reference has a non-zero H matrix element
     * @return
     */
    Walker& create_row_(uint_t icycle, const field::Mbf &mbf, ham_comp_t hdiag, const v_t<bool>& refconns);


    Walker& create_row_(uint_t icycle, const field::Mbf &mbf,
                        const ham_comp_t &hdiag, bool refconn) {
        return create_row_(icycle, mbf, hdiag, v_t<bool>(npart(), refconn));
    }

    /**
     * Called on all ranks, dispatching create_row_ on the assigned rank only
     */
    TableBase::Loc create_row(uint_t icycle, const field::Mbf &mbf, ham_comp_t hdiag, const v_t<bool>& refconns);


    TableBase::Loc create_row(uint_t icycle, const field::Mbf &mbf, ham_comp_t hdiag, bool refconn) {
        return create_row(icycle, mbf, hdiag, v_t<bool>(npart(), refconn));
    }

    Spawn& add_spawn(const field::Mbf &dst_mbf, wf_t delta, bool initiator, bool deterministic, uint_t dst_ipart);

    Spawn& add_spawn(const field::Mbf &dst_mbf, wf_t delta, bool initiator, bool deterministic, uint_t dst_ipart,
                     const field::Mbf &src_mbf, wf_t src_weight);

    uint_t npart() const {
        return m_format.m_nelement;
    }

    void fci_init(const Hamiltonian& h, FciInitOptions opts, uint_t max_ncomm=1000ul);
private:


    void orthogonalize(NdReduction<wf_t, 3>& overlaps,
                       uint_t iroot, uint_t jroot, uint_t ireplica) {
        ASSERT(iroot<=jroot);
        auto& row = m_store.m_row;
        const auto ipart_src = m_format.flatten({iroot, ireplica});
        const auto ipart_dst = m_format.flatten({jroot, ireplica});
        overlaps.m_local[{iroot, jroot, ireplica}] +=
                arith::conj(row.m_weight[ipart_src]) * row.m_weight[ipart_dst];
        if (jroot+1<nroot()) {
            // there is another part to project onto
            const auto ipart_next = m_format.flatten({jroot+1, ireplica});
            overlaps.m_local[{iroot, jroot + 1, ireplica}] +=
                    arith::conj(row.m_weight[ipart_src]) * row.m_weight[ipart_next];
        }
        if (iroot<jroot){
            const auto& overlap = overlaps.m_reduced[{iroot, jroot, ireplica}];
            const auto& norm = overlaps.m_reduced[{iroot, iroot, ireplica}];
            ASSERT(std::abs(norm)>1e-12);
            const auto gs_coeff = overlap / norm;
            change_weight(row, ipart_dst, -gs_coeff*row.m_weight[ipart_src]);
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
        NdReduction<wf_t, 3> overlaps({nroot(), nroot(), nreplica()});
        auto& row = m_store.m_row;
        for (uint_t iroot = 0ul; iroot < nroot(); ++iroot) {
            for (uint_t jroot = iroot; jroot < nroot(); ++jroot) {
                for (uint_t ireplica = 0ul; ireplica < nreplica(); ++ireplica) {
                    for (row.restart(); row; ++row) {
                        if (!row.m_mbf.is_zero()) orthogonalize(overlaps, iroot, jroot, ireplica);
                    }
                    overlaps.all_sum();
                }
            }
        }
    }

    void save(const hdf5::NodeWriter &parent) const {
        auto& row = m_store.m_row;
        hdf5::GroupWriter gw(parent, "wf");
        row.m_mbf.save(gw, true);
        row.m_weight.save(gw, true);
    }
};

#endif //M7_WAVEFUNCTION_H
