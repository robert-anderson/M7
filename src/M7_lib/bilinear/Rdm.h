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

class Rdm : public communicator::MappedSend<MaeRow, MaeRow> {
public:
    const sys::Size m_basis_size;
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

protected:
    /**
     * enumerators of the promotions of contributing excitation signatures to the ranksig of the RDM
     */
    v_t<FermionPromoter> m_frm_promoters;
    /**
     * indices of the full second quantised string
     */
    buffered::MaeInds m_full_inds;
    /**
     * indices of any contracted intermediate, not the full second quantised string
     */
    buffered::MaeInds m_uncontracted_inds;
    const str_t m_name;

    /**
     * objects for table traversal
     */
    MaeRow m_send_row;
    MaeRow m_recv_row;
    MaeRow m_store_row;

    static uint_t nrow_estimate(uint_t nfrm_cre, uint_t nfrm_ann, uint_t nbos_cre, uint_t nbos_ann, sys::Size basis_size);

    static uint_t nrow_estimate(uint_t exsig, sys::Size basis_size);

    str_t name(str_t str, uint_t ranksig) const;

    void add_to_send_table(const field::MaeInds& inds, wf_t contrib);

    virtual void frm_make_contribs(const field::FrmOnv& src_onv, const conn::FrmOnv& conn,
                                   const com_ops::Frm& com, const wf_t& contrib);

    virtual void frmbos_make_contribs(const field::FrmBosOnv& src_onv, const conn::FrmBosOnv& conn,
                                      const com_ops::FrmBos& com, const wf_t& contrib);

    virtual void bos_make_contribs(const field::BosOnv& /*src_onv*/, const conn::BosOnv& /*conn*/,
                                   const com_ops::Bos& /*com*/, const wf_t& /*contrib*/) {
        ABORT("not yet implemented");
    }

    static uint_t nrec_est(sys::Size basis_size, uint_t indsig) {
        const auto nrec_frm = integer::combinatorial(basis_size.m_frm.m_nspinorb, decode_nfrm_cre(indsig));
        return nrec_frm * decode_nbos(indsig);
    }

public:

    str_t name() const {
        return name(m_name, m_ranksig);
    }

    const uint_t m_nelec;

    Rdm(uint_t ranksig, uint_t indsig, sys::Size basis_size, uint_t nelec, uint_t nvalue,
        DistribOptions dist_opts, Sizing store_sizing, Sizing comm_sizing, str_t name="");

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
    Rdm(const conf::Rdms& opts, uint_t ranksig, uint_t indsig,
        sys::Size basis_size, uint_t nelec, uint_t nvalue, str_t name=""):
        Rdm(ranksig, indsig, basis_size, nelec, nvalue, opts.m_distribution,
            // store sizing
            Sizing{nrec_est(basis_size, indsig), opts.m_buffers.m_store_exp_fac},
            // send/recv sizing
            Sizing{nrec_est(basis_size, indsig), opts.m_buffers.m_comm_exp_fac}, name){}

    void end_cycle();

    void save(hdf5::NodeWriter& gw) const;

    void make_contribs(const field::FrmOnv& src_onv, const conn::FrmOnv& conn,
                       const com_ops::Frm& com, const wf_t& contrib) {
        frm_make_contribs(src_onv, conn, com, contrib);
    }

    void make_contribs(const field::FrmBosOnv& src_onv, const conn::FrmBosOnv& conn,
                       const com_ops::FrmBos& com, const wf_t& contrib) {
        frmbos_make_contribs(src_onv, conn, com, contrib);
    }

    void make_contribs(const field::BosOnv& src_onv, const conn::BosOnv& conn,
                       const com_ops::Bos& com, const wf_t& contrib) {
        bos_make_contribs(src_onv, conn, com, contrib);
    }
};

class FockRdm4 : public Rdm {

    dense::SquareMatrix<ham_t> m_fock;

public:
    FockRdm4(const conf::Rdms& opts, sys::Size basis_size, uint_t nelec, uint_t nvalue);

    /**
     * override the default method to implement on-the-fly contraction
     */
    void frm_make_contribs(const field::FrmOnv& src_onv, const conn::FrmOnv& conn,
                           const FrmOps& com, const wf_t& contrib) override;
};

class Rdms : public Archivable {
    /**
     * RDM objects managed by this instance
     */
    std::array<std::unique_ptr<Rdm>, exsig::c_ndistinct> m_rdms;
    /**
     * if true, the spinfree versions of the full RDMs are all computed and output to the HDF5 archive (note that this
     * is only a finalization procedure, NOT done on the fly)
     */
    const bool m_spinfree;
    /**
     * optionally-allocatable contraction of the 4RDM with the active space generalize Fock matrix for use in CASPT2
     */
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

    void save_fn(const hdf5::NodeWriter& parent) override;
};

#endif //M7_RDM_H
