//
// Created by Robert J. Anderson on 11/08/2021.
//

#ifndef M7_RDM_H
#define M7_RDM_H

#include <M7_lib/basis/Suites.h>
#include <M7_lib/wavefunction/Reference.h>
#include <M7_lib/propagator/Propagator.h>
#include <M7_lib/mae/MaeTable.h>
#include <M7_lib/field/Fields.h>
#include <M7_lib/io/Archivable.h>
#include <M7_lib/hdf5/Group.h>

#include "FermionPromoter.h"

using namespace exsig;

class Rdm : public Communicator<MaeRow, MaeRow, true> {
protected:
    /**
     * rank signature of the RDM, along with convenient decoded
     */
    const uint_t m_ranksig;
    const uint_t m_rank, m_nfrm_cre, m_nfrm_ann, m_nbos_cre, m_nbos_ann;
    /**
     * index signature (different to ranksig if RDM is stored as an on-the-fly contraction)
     */
    const uint_t m_indsig;
    const uint_t m_rank_ind, m_nfrm_cre_ind, m_nfrm_ann_ind, m_nbos_cre_ind, m_nbos_ann_ind;
    /**
     * enumerators of the promotions of contributing excitation signatures to the ranksig of the RDM
     */
    v_t<FermionPromoter> m_frm_promoters;
    /**
     * indices of the full second quantised string
     */
    buffered::MaeInds m_full_inds;
    const str_t m_name;

    static uint_t nrow_estimate(uint_t nfrm_cre, uint_t nfrm_ann, uint_t nbos_cre, uint_t nbos_ann, sys::Size basis_size);

    static uint_t nrow_estimate(uint_t exsig, sys::Size basis_size);

    str_t name(str_t str, uint_t ranksig) const {
        return str.empty() ? exsig::to_string(ranksig) : str;
    }

public:

    str_t name() const {
        return name(m_name, m_ranksig);
    }

    const uint_t m_nelec;

    /**
     * @param opts
     *  RDM section of the config document
     * @param ranksig
     *  rank of the SQ operators in each contribution
     * @param basis_size
     *  dimensions of the stored basis
     * @param nelec
     *  number of electrons to use in enforcing probability-conserving trace
     *  TODO: generalize to use sys::Particles
     * @param nvalue
     *  number of values to encode in each RDM element
     * @param name
     *  string identifier for logging and archive output. if empty, this is generated from the ranksig
     * @param indsig
     *  number of each species of SQ operator to store in the structure (equal to ranksig for ordinary, uncontracted RDMs)
     */
    Rdm(const conf::Rdms& opts, uint_t ranksig, uint_t indsig, sys::Size basis_size, uint_t nelec, uint_t nvalue, str_t name="");

    void make_contribs(const field::FrmOnv& src_onv, const conn::FrmOnv& conn,
                       const com_ops::Frm& com, const wf_t& contrib);

    void make_contribs(const field::FrmBosOnv& src_onv, const conn::FrmBosOnv& conn,
                       const com_ops::FrmBos& com, const wf_t& contrib);

    void make_contribs(const field::BosOnv& /*src_onv*/, const conn::BosOnv& /*conn*/,
                       const com_ops::Bos& /*com*/, const wf_t& /*contrib*/) {
        ABORT("not yet implemented");
    }

    void end_cycle();

    void save(hdf5::NodeWriter& gw) const;
};

class FockRdm4 : public Rdm {

    /**
     * indices of the intermediate, not the full second quantised string
     */
    buffered::MaeInds m_uncontracted_inds;
    dense::SquareMatrix<ham_t> m_fock;

public:
    FockRdm4(const conf::Rdms& opts, sys::Size basis_size, uint_t nelec, uint_t nvalue):
        Rdm(opts, exsig::ex_4400, exsig::ex_3300, basis_size, nelec, nvalue, "4400f"),
        m_uncontracted_inds(ex_3300), m_fock(basis_size.m_frm.m_nsite){
        logging::info("loading generalized Fock matrix for CASPT2 contracted 4RDM accumulation");
        //opts.m_fock_4rdm.m_fock_path;
        hdf5::FileReader reader(opts.m_fock_4rdm.m_fock_path);
        REQUIRE_TRUE(reader.child_exists("ACT_FOCK_INDEX"), "invalid fock matrix file contents");
        REQUIRE_TRUE(reader.child_exists("ACT_FOCK_VALUES"), "invalid fock matrix file contents");
        v_t<int64_t> inds;
        reader.read_data("ACT_FOCK_INDEX", inds);
        v_t<double> values;
        reader.read_data("ACT_FOCK_VALUES", values);
        REQUIRE_EQ(inds.size(), values.size()*2, "incorrect number of indices");
        for (uint_t i = 0ul; i<values.size(); ++i) {
            m_fock(inds[2*i]-1, inds[2*i+1]-1) = values[i];
        }
    }

    void make_contribs(const field::FrmOnv& src_onv, const conn::FrmOnv& conn, const FrmOps& com, const wf_t& contrib) {
        // TODO: use virtual method to abstract the promotion process away from the individual contributions
        const auto exlvl = conn.m_cre.size();
        DEBUG_ASSERT_TRUE(conn.m_ann.size() <= m_nfrm_ann && conn.m_cre.size() <= m_nfrm_cre,
                          "this method should not have been delegated given the exsig of the contribution");
        /*
         * number of "inserted" fermion creation/annihilation operator pairs
         */
        const auto nins = m_rank - exlvl;
        /*
         * this determines the precomputed promoter required
         */
        const auto& promoter = m_frm_promoters[nins];
        /*
         * apply each combination of the promoter deterministically
         */
        for (uint_t icomb = 0ul; icomb < promoter.m_ncomb; ++icomb) {
            auto phase = promoter.apply(icomb, conn, com, m_full_inds.m_frm);
            /*
             * include the Fermi phase of the excitation
             */
            phase = phase ^ conn.phase(src_onv);
            bool contract_phase = true;
            // TODO: diagonal Fock optimisation
            for (uint_t icre_contract = 0ul; icre_contract < 4ul; ++icre_contract){
                const auto isite_cre = src_onv.m_basis.isite(m_full_inds.m_frm.m_cre[icre_contract]);
                for (uint_t iann_contract = 0ul; iann_contract < 4ul; ++iann_contract) {
                    const auto isite_ann = src_onv.m_basis.isite(m_full_inds.m_frm.m_ann[iann_contract]);

                    const auto fock_element = m_fock(isite_cre, isite_ann);

                    contract_phase = !contract_phase;
                    /*
                     * fill the uncontracted indices (those identifying the elements of the intermediate)
                     */
                    for (uint_t iuncontract = 0ul; iuncontract < 4ul; ++iuncontract) {
                        uint_t i;
                        if (iuncontract!=icre_contract) {
                            i = iuncontract - (iuncontract > icre_contract ? 1 : 0);
                            m_uncontracted_inds.m_frm.m_cre[i] = m_full_inds.m_frm.m_cre[iuncontract];
                        }
                        if (iuncontract!=iann_contract) {
                            i = iuncontract - (iuncontract > iann_contract ? 1 : 0);
                            m_uncontracted_inds.m_frm.m_ann[i] = m_full_inds.m_frm.m_ann[iuncontract];
                        }
                    }

                    auto irank_send = irank(m_uncontracted_inds);
                    DEBUG_ASSERT_TRUE(m_full_inds.is_ordered(),
                                      "operators of each kind should be stored in ascending order of their orbital (or mode) index");
                    auto& send_table = send(irank_send);
                    uint_t irow = *send_table[m_uncontracted_inds];
                    if (irow == ~0ul) irow = send_table.insert(m_uncontracted_inds);
                    send_table.m_row.jump(irow);
                    /*
                     * include the Fermi phase of the rearrangement
                     */
                    send_table.m_row.m_values[0] += ((phase ^ contract_phase) ? -contrib : contrib)*fock_element;
                }
            }
        }
    }
};

class Rdms : public Archivable {
    std::array<std::unique_ptr<Rdm>, exsig::c_ndistinct> m_rdms;
    std::unique_ptr<FockRdm4> m_fock_rdm4;
    const uintv_t m_rdm_ranksigs;
    const std::array<uintv_t, exsig::c_ndistinct> m_exsig_ranks;

    suite::Conns m_work_conns;
    suite::ComOps m_work_com_ops;

    std::array<uintv_t, exsig::c_ndistinct> make_exsig_ranks() const;

public:
    const bool m_explicit_ref_conns;
    const Epoch& m_accum_epoch;
    Reduction<wf_t> m_total_norm;
    const uint_t m_nelec;

    Rdms(const conf::Rdms& opts, uintv_t ranksigs, sys::Size basis_size, uint_t nelec, const Epoch& accum_epoch);

    operator bool() const;

    bool takes_contribs_from(uint_t exsig) const;

    void make_contribs(const field::Mbf& src_onv, const conn::Mbf& conn,
                       const com_ops::Mbf& com, const wf_t& contrib);

    void make_contribs(const field::Mbf& src_onv, const field::Mbf& dst_onv, const wf_t& contrib);

    void make_contribs(const SpawnTableRow& recv_row, const WalkerTableRow& dst_row, const Propagator& prop);

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
     * By the time MAE contributions are being made, the death step has already been applied, and so the pre-
     * death value of the weight must be reconstituted by undoing the scaling, thus, the pre-death value of Cdst
     * is just Cdst/(1 - tau (Hii-shift))
     *
     * @param src_mbf
     * @param dst_row
     * @param prop
     * @param refs
     * @param ipart_dst
     */
//    void make_contribs(const field::Mbf& src_mbf, const wf_t& src_weight, const WalkerTableRow& dst_row,
//                       const Propagator& prop, const References& refs, const uint_t& ipart_dst) {
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

    /**
     * @param ham
     *  hamiltonian corresponding to the evolution of the wavefunction(s) used in the RDM estimation
     * @return
     *  true only if the ranks of RDMs estimated are sufficient for pseudo-variational energy estimation via contraction
     *  of the RDMs with the Hamiltonian coefficients.
     */
    bool is_energy_sufficient(const Hamiltonian& ham) const;

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
    ham_comp_t get_energy(const FrmHam& ham) const;

    /**
     * compute the RDM energy contribution from the boson number-nonconserving terms
     * @param ham
     *  boson ladder-operator (pure and coupled) hamiltonian
     * @return
     */
    ham_comp_t get_energy(const FrmBosHam& ham, uint_t nelec, uint_t exsig) const;

    ham_comp_t get_energy(const FrmBosHam& ham, uint_t nelec) const {
        return get_energy(ham, nelec, exsig::ex_1101) + get_energy(ham, nelec, exsig::ex_1110);
    }

    /**
     * compute the RDM energy contribution from the boson number-conserving terms
     * @param ham
     *  boson number conserving hamiltonian
     * @return
     */
    ham_comp_t get_energy(const BosHam& ham) const;

    /**
     * @param ham
     *  full hamiltonian including all terms
     * @return
     *  E_RDM = E_2RDM + E_RDM_ladder + E_RDM_boson
     */
    ham_comp_t get_energy(const Hamiltonian& ham) const {
        if (!is_energy_sufficient(ham)) return 0.0;
        return get_energy(ham.m_frm) + get_energy(ham.m_frmbos, m_nelec) + get_energy(ham.m_bos);
    }

private:
    void load_fn(const hdf5::NodeReader& /*parent*/) override {

    }

    void save_fn(const hdf5::NodeWriter& parent) override {
        if (!m_accum_epoch) {
            logging::warn("MAE accumulation epoch was not reached in this calculation: omitting RDM save");
            return;
        }
        hdf5::GroupWriter gw(parent, "rdms");
        gw.write_data("norm", m_total_norm.m_reduced);
        for (const auto& i: m_rdm_ranksigs) {
            DEBUG_ASSERT_TRUE(m_rdms[i].get(), "active ranksig was not allocated!");
            m_rdms[i]->save(gw);
        }
        if (m_fock_rdm4) m_fock_rdm4->save(gw);
    }
};

#endif //M7_RDM_H
