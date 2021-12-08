//
// Created by rja on 11/08/2021.
//

#ifndef M7_RDM_H
#define M7_RDM_H

#include <src/core/basis/Suites.h>
#include <src/core/wavefunction/Reference.h>
#include <src/core/dynamics/Propagator.h>
#include "src/core/mae/MaeTable.h"
#include "src/core/field/Fields.h"
#include "src/core/table/Communicator.h"
#include "src/core/io/Archivable.h"
#include "FermionPromoter.h"

using namespace exsig_utils;

class Rdm : public Communicator<MaeRow, MaeRow, true> {
    const size_t m_ranksig;
    const size_t m_rank, m_nfrm_cre, m_nfrm_ann, m_nbos_cre, m_nbos_ann;
    std::vector<FermionPromoter> m_frm_promoters;
    buffered::MaeInds m_lookup_inds;
    static size_t nrow_estimate(size_t nfrm_cre, size_t nfrm_ann, size_t nbos_cre, size_t nbos_ann, BasisDims bd);

    static size_t nrow_estimate(size_t exsig, BasisDims bd);

public:
    const BasisDims m_bd;

    Rdm(const fciqmc_config::Rdms &opts, size_t ranksig, BasisDims bd, size_t nelec, size_t nvalue);

    void make_contribs(const field::FrmOnv &src_onv, const conn::FrmOnv &conn,
                       const com_ops::Frm &com, const defs::wf_t &contrib);

    void make_contribs(const field::FrmBosOnv &src_onv, const conn::FrmBosOnv &conn,
                       const com_ops::FrmBos &com, const defs::wf_t &contrib);

    void make_contribs(const field::BosOnv &src_onv, const conn::BosOnv &conn,
                       const com_ops::Bos &com, const defs::wf_t &contrib) {
        ABORT("not yet implemented");
    }

    void end_cycle();

    void save(hdf5::GroupWriter &gw) const;
};

class Rdms : public Archivable {
    std::array<std::unique_ptr<Rdm>, defs::nexsig> m_rdms;
    const defs::inds m_active_ranksigs;
    const std::array<defs::inds, defs::nexsig> m_exsig_ranks;

    suite::Conns m_work_conns;
    suite::ComOps m_work_com_ops;

    std::array<defs::inds, defs::nexsig> make_exsig_ranks() const;

public:
    const bool m_explicit_ref_conns;
    const Epoch &m_accum_epoch;
    Reduction<defs::wf_comp_t> m_total_norm;

    Rdms(const fciqmc_config::Rdms &opts, defs::inds ranksigs, BasisDims bd, size_t nelec, const Epoch &accum_epoch);

    operator bool() const;

    bool takes_contribs_from(const size_t &exsig) const;

    void make_contribs(const field::Mbf &src_onv, const conn::Mbf &conn,
                       const com_ops::Mbf &com, const defs::wf_t &contrib);

    void make_contribs(const field::Mbf &src_onv, const field::Mbf &dst_onv, const defs::wf_t &contrib);

    void make_contribs(const SpawnTableRow &recv_row, const WalkerTableRow &dst_row, const Propagator &prop);

    /**
     * TODO: remove. this is now done in the Annihilator class
     * We need to be careful of the intermediate state of the walker weights.
     * if src_weight is taken from the wavefunction at cycle i, dst_weight is at an intermediate value equal to
     * the wavefunction at cycle i with the diagonal part of the propagator already applied. We don't want to
     * introduce a second post-annihilation loop over occupied MBFs to apply the diagonal part of the
     * propagator, so for MEVs, the solution is to reconstitute the value of the walker weight before the
     * diagonal death/cloning.
     *
     * The death-step behaviour of the exact (and stochastic on average) propagator is to scale the WF:
     * Ci -> Ci*(1 - tau (Hii-shift)).
     * By the time MEV contributions are being made, the death step has already been applied, and so the pre-
     * death value of the weight must be reconstituted by undoing the scaling, thus, the pre-death value of Cdst
     * is just Cdst/(1 - tau (Hii-shift))
     *
     * @param src_mbf
     * @param dst_row
     * @param prop
     * @param refs
     * @param ipart_dst
     */
//    void make_contribs(const field::Mbf &src_mbf, const defs::wf_t &src_weight, const WalkerTableRow &dst_row,
//                       const Propagator &prop, const References &refs, const size_t &ipart_dst) {
//        if (!*this) return;
//        if (!m_accum_epoch) return;
//
//        auto ipart_dst_replica = dst_row.ipart_replica(ipart_dst);
//        double dupl_fac = 1.0 / dst_row.nreplica();
//
//        auto dst_weight_before_death = dst_row.m_weight[ipart_dst_replica];
//        dst_weight_before_death /= 1.0 - prop.tau() * (dst_row.m_hdiag - prop.m_shift.m_values[ipart_dst_replica]);
//        make_contribs(src_mbf, dst_row.m_mbf, dupl_fac * src_weight * dst_weight_before_death);
//    }


    bool all_stores_empty() const;

    void end_cycle();

    bool is_energy_sufficient(const Hamiltonian &ham) const;

    /**
     * compute the fermion 2-RDM energy
     * @param ham
     * @return
     *  MPI-reduced sum of 0, 1, and 2 body parts of the RDM energy
     *
     *  E_RDM = h0 + h1[i,j] * rdm1[i,j] + <ij|kl> * rdm2[i,j,k,l]
     *
     *  rdm1[i,j] = sum_k rdm2[i,k,j,k] / (n_elec - 1)
     */
    defs::ham_comp_t get_energy(const FermionHamiltonian &ham) const;

    /**
     * compute the RDM energy contribution from the boson number-nonconserving terms
     * @param ham
     *  boson ladder-operator (pure and coupled) hamiltonian
     * @return
     */
    defs::ham_comp_t get_energy(const LadderHamiltonian &ham, size_t nelec, size_t exsig) const;

    defs::ham_comp_t get_energy(const LadderHamiltonian &ham, size_t nelec) const {
        return get_energy(ham, nelec, exsig_utils::ex_1101) + get_energy(ham, nelec, exsig_utils::ex_1110);
    }

    /**
     * compute the RDM energy contribution from the boson number-conserving terms
     * @param ham
     *  boson number conserving hamiltonian
     * @return
     */
    defs::ham_comp_t get_energy(const BosonHamiltonian &ham) const;

    /**
     * @param ham
     *  full hamiltonian including all terms
     * @return
     *  E_RDM = E_2RDM + E_RDM_ladder + E_RDM_boson
     */
    defs::ham_comp_t get_energy(const Hamiltonian &ham) const {
        if (!is_energy_sufficient(ham)) return 0.0;
        return get_energy(*ham.m_frm) + get_energy(ham.m_ladder, ham.nelec()) + get_energy(ham.m_bos);
    }

private:
    void load_fn(hdf5::GroupReader &parent) override {

    }

    void save_fn(hdf5::GroupWriter &parent) override {
        if (!m_accum_epoch) {
            log::warn("MAE accumulation epoch was not reached in this calculation: omitting RDM save");
            return;
        }
        hdf5::GroupWriter gw("rdms", parent);
        gw.save("norm", m_total_norm.m_reduced);
        for (const auto &i: m_active_ranksigs) {
            DEBUG_ASSERT_TRUE(m_rdms[i].get(), "active ranksig was not allocated!");
            m_rdms[i]->save(gw);
        }
    }
};

#endif //M7_RDM_H
