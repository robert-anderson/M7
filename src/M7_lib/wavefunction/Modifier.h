//
// Created by rja on 24/02/23.
//

#ifndef M7_MODIFIER_H
#define M7_MODIFIER_H

#include "Reference.h"

#if 0
namespace wf {

    /**
     * when making changes to the FCI wavefunction, we must:
     *  1. log changes to the walker weights, number of initiators, and other statistics
     *  2. keep up-to-date the quantities computed and cached within the Walker rows of m_wf.m_store
     * the purpose of this abstraction is to provide a single object to:
     *  1. statefully log these changes
     *  2. allow cached quantities to be updated without having the wavefunction class make explicit reference to all
     *     objects necessary to compute these updates. Most importantly: we prefer that the Wavefunction does not
     *     depend on the Hamiltonian object, even though it caches the diagonal H elements of all stored walkers. This
     *     is neatly achieved by implementing the following modification helper
     */
    struct Modifier {
        /**
         * solution vector storing multiple eigenvectors of H in distributed memory
         */
        wf::Fci &m_wf;
        /**
         * reference many-body basis functions (MBFs)
         */
        wf::Refs &m_refs;
        /**
         * Hamiltonian under which m_wf is being evolved by some Propagator
         */
        const Hamiltonian &m_ham;

        /**
         * collection of all reductions which are summed at the end of every cycle
         */
        reduction::Syndicate m_summables;

        /**
         * number of initiator MBFs in each part of the WF
         */
        reduction::NdArray<uint_t, c_ndim_wf> m_ninitiator;
        /**
         * number of initiator MBFs in each part of the WF due to permanitiator status
         */
        reduction::NdArray<uint_t, c_ndim_wf> m_ninitiator_perma;
        /**
         * number of MBFs with any associated weight in any part
         */
        reduction::Scalar<uint_t> m_nocc_mbf;
        /**
         * change in the number of occupied MBFs
         */
        reduction::Scalar<int> m_delta_nocc_mbf;
        /**
         * L1 norm of each part of the WF
         */
        reduction::cyclic::NdArray<wf_comp_t, c_ndim_wf> m_nwalker;
        /**
         * square of the L2 norm of each part of the WF
         */
        reduction::NdArray<wf_comp_t, c_ndim_wf> m_l2_norm_square;
        /**
         * change in the L2 norm
         */
        reduction::NdArray<wf_comp_t, c_ndim_wf> m_delta_l2_norm_square;
        /**
         * number of walkers received in spawning process
         */
        reduction::NdArray<wf_comp_t, c_ndim_wf> m_nspawned;
        /**
         * number of walkers annihilated in the loop_over_spawned method for each part
         */
        reduction::NdArray<wf_comp_t, c_ndim_wf> m_nannihilated;

        Modifier(wf::Fci& wf, wf::Refs& refs, const Hamiltonian& ham):
            m_wf(wf), m_refs(refs), m_ham(ham),
            m_ninitiator(m_wf.m_format),
            m_ninitiator_perma(m_wf.m_format),
            m_nwalker(m_wf.m_format),
            m_l2_norm_square(m_wf.m_format),
            m_delta_l2_norm_square(m_wf.m_format),
            m_nspawned(m_wf.m_format),
            m_nannihilated(m_wf.m_format) {
            m_summables.add_members(m_ninitiator, m_nocc_mbf, m_delta_nocc_mbf,
                                    m_nwalker, m_l2_norm_square, m_delta_l2_norm_square,
                                    m_nspawned, m_nannihilated);
        }

        void set_weight(Walker& walker, uint_t ipart, wf_t new_weight) {
            wf_t& weight = walker.m_weight[ipart];
            m_nwalker.
            m_delta_nwalker.m_local[ipart] += std::abs(new_weight);
            m_delta_nwalker.m_local[ipart] -= std::abs(weight);
            m_delta_l2_norm_square.m_local[ipart] += std::pow(std::abs(new_weight), 2.0);
            m_delta_l2_norm_square.m_local[ipart] -= std::pow(std::abs(weight), 2.0);
            weight = new_weight;
        }

        void change_weight(Walker& walker, uint_t ipart, wf_t delta) {
            set_weight(walker, ipart, walker.m_weight[ipart] + delta);
        }

        void scale_weight(Walker& walker, uint_t ipart, double factor) {
            set_weight(walker, ipart, factor * walker.m_weight[ipart]);
        }

        void zero_weight(Walker& walker, uint_t ipart) {
            set_weight(walker, ipart, 0.0);
        }

        void remove_row(Walker& walker) {
            DEBUG_ASSERT_TRUE(m_store.lookup(walker.m_mbf), "MBF doesn't exist in table!");
            for (uint_t ipart = 0ul; ipart < m_wf.m_format.m_nelement; ++ipart) {
                zero_weight(walker, ipart);
                m_delta_nocc_mbf.m_local--;
            }
            m_wf.m_store.erase(walker.m_mbf);
        }

        Walker& create_row_(uint_t icycle, const Mbf& mbf, ham_comp_t hdiag, const v_t<bool>& refconns) {
            DEBUG_ASSERT_EQ(refconns.size(), npart(), "should have as many reference rows as WF parts");
            DEBUG_ASSERT_TRUE(mpi::i_am(m_dist.irank(mbf)),
                              "this method should only be called on the rank responsible for storing the MBF");
            auto& row = m_wf.m_store.insert(mbf);
            m_delta_nocc_mbf.m_local++;
            DEBUG_ASSERT_EQ(row.key_field(), mbf, "MBF was not properly copied into key field of WF row");
            row.m_hdiag = hdiag;
            for (uint_t ipart=0ul; ipart < npart(); ++ipart)
                row.m_ref_conn.put(ipart, refconns[ipart]);
            /*
             * we need to be very careful here of off-by-one-like mistakes. the initial walker is "created" at the beginning
             * of MC cycle 0, and so the stats line output for cycle 0 will show that the number of walkers is the initial
             * occupation of the initial row. if a spawning event leads to the creation of another row, it is created on
             * iteration 1 even though it is added in the annihilating call of iteration 0. so, if this method is called in
             * the annihilating process of MC cycle i, it actually "becomes occupied" on cycle i+1.
             */
            if (storing_av_weights()) {
                row.m_icycle_occ = icycle+1;
                row.m_average_weight = 0;
            }
            return row;
        }

        TableBase::Loc wf::Fci::create_row(uint_t icycle, const Mbf& mbf, ham_comp_t hdiag, const v_t<bool>& refconns) {
            const uint_t irank = m_dist.irank(mbf);
            uint_t irec;
            if (mpi::i_am(irank)) {
                irec = create_row_(icycle, mbf, hdiag, refconns).index();
            }
            mpi::bcast(irec, irank);
            return {irank, irec};
        }

        Spawn& wf::Fci::add_spawn(const field::Mbf& dst_mbf, wf_t delta, bool initiator,
                                  bool deterministic, uint_t dst_ipart) {
            auto& dst_table = send(m_dist.irank(dst_mbf));

            auto& spawn = dst_table.m_row;
            spawn.push_back_jump();

            spawn.m_dst_mbf = dst_mbf;
            spawn.m_delta_weight = delta;
            spawn.m_src_initiator = initiator;
            spawn.m_src_deterministic = deterministic;
            spawn.m_ipart_dst = dst_ipart;
            return spawn;
        }

        Spawn& wf::Fci::add_spawn(const field::Mbf& dst_mbf, wf_t delta, bool initiator, bool deterministic,
                                  uint_t dst_ipart, const field::Mbf& src_mbf, wf_t src_weight) {
            auto& spawn = add_spawn(dst_mbf, delta, initiator, deterministic, dst_ipart);
            if (spawn.m_send_parents) {
                spawn.m_src_mbf = src_mbf;
                spawn.m_src_weight = src_weight;
            }
            DEBUG_ASSERT_NE(dst_mbf, src_mbf, "spawning diagonally");
            return spawn;
        }


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
            for (uint_t ipart = 0ul; ipart < m_format.m_nelement; ++ipart) set_weight(walker, ipart, new_weight);
        }

        void set_weight(Walker& walker, const Numbers<wf_t, c_ndim_wf>& new_weight) {
            for (uint_t i = 0ul; i < m_format.m_nelement; ++i) set_weight(walker, i, new_weight[i]);
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
        Walker& create_row_(uint_t icycle, const Mbf& mbf, ham_comp_t hdiag, const v_t<bool>& refconns);


        Walker& create_row_(uint_t icycle, const Mbf& mbf,
                            const ham_comp_t& hdiag, bool refconn) {
            return create_row_(icycle, mbf, hdiag, v_t<bool>(npart(), refconn));
        }

        /**
         * Called on all ranks, dispatching create_row_ on the assigned rank only
         */
        TableBase::Loc create_row(uint_t icycle, const Mbf& mbf, ham_comp_t hdiag, const v_t<bool>& refconns);


        TableBase::Loc create_row(uint_t icycle, const Mbf& mbf, ham_comp_t hdiag, bool refconn) {
            return create_row(icycle, mbf, hdiag, v_t<bool>(npart(), refconn));
        }

        Spawn& add_spawn(const Mbf& dst_mbf, wf_t delta, bool initiator, bool deterministic, uint_t dst_ipart);

        Spawn& add_spawn(const Mbf& dst_mbf, wf_t delta, bool initiator, bool deterministic, uint_t dst_ipart,
                         const Mbf& src_mbf, wf_t src_weight);

    private:


        void orthogonalize(NdReduction<wf_t, 3>& overlaps,
                           uint_t iroot, uint_t jroot, uint_t ireplica) {
            ASSERT(iroot <= jroot);
            auto& row = m_store.m_row;
            const auto ipart_src = m_format.flatten({iroot, ireplica});
            const auto ipart_dst = m_format.flatten({jroot, ireplica});
            overlaps.m_local[{iroot, jroot, ireplica}] +=
                    arith::conj(row.m_weight[ipart_src]) * row.m_weight[ipart_dst];
            if (jroot + 1 < nroot()) {
                // there is another part to project onto
                const auto ipart_next = m_format.flatten({jroot + 1, ireplica});
                overlaps.m_local[{iroot, jroot + 1, ireplica}] +=
                        arith::conj(row.m_weight[ipart_src]) * row.m_weight[ipart_next];
            }
            if (iroot < jroot) {
                const auto& overlap = overlaps.m_reduced[{iroot, jroot, ireplica}];
                const auto& norm = overlaps.m_reduced[{iroot, iroot, ireplica}];
                ASSERT(std::abs(norm) > 1e-12);
                const auto gs_coeff = overlap / norm;
                change_weight(row, ipart_dst, -gs_coeff * row.m_weight[ipart_src]);
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
    };
}

#endif //M7_MODIFIER_H
#endif //M7_MODIFIER_H
