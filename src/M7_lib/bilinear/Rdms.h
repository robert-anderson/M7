//
// Created by rja on 03/11/22.
//

#ifndef M7_RDMS_H
#define M7_RDMS_H

#include "Rdm.h"
#include "FockRdm4.h"

class Rdms : public Archivable {
    /**
     * RDM objects managed by this instance
     */
    std::forward_list<std::unique_ptr<Rdm>> m_rdms;
    /**
     * if true, the spinfree versions of the full RDMs are all computed and output to the HDF5 archive (note that this
     * is only a finalization procedure, NOT done on the fly)
     */
    const bool m_spinfree;
    /**
     * each array element is indexed by an exsig, and contains a list of pointers to RDM objects which take
     * contributions from that exsig, either because the RDM's ranksig matches the exsig, or the exsig is promotable to
     * it
     */
    typedef std::array<std::forward_list<Rdm*>, exsig::c_ndistinct> exsig_to_rdms_t;
    exsig_to_rdms_t m_exsig_to_rdms;

    /**
     * each array element is indexed by a ranksig, and points to the pure RDM of that rank
     */
    typedef std::array<Rdm*, exsig::c_ndistinct> pure_rdms_t;
    pure_rdms_t m_pure_rdms;

    suite::Conns m_work_conns;
    suite::ComOps m_work_com_ops;

    exsig_to_rdms_t make_exsig_to_rdms() const;

public:
    const bool m_explicit_ref_conns;
    const Epoch& m_accum_epoch;
    Reduction<wf_t> m_total_norm;

    Rdms(const conf::Rdms& opts, uintv_t ranksigs, sys::Sector sector, const Epoch& accum_epoch);

    operator bool() const;

    bool takes_contribs_from(uint_t exsig) const;

    void make_contribs(const field::Mbf& src_onv, const conn::Mbf& conn,
                       const com_ops::Mbf& com, const wf_t& contrib);

    void make_contribs(const field::Mbf& src_onv, const field::Mbf& dst_onv, const wf_t& contrib);

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
    ham_comp_t get_energy(const FrmBosHam& ham, uint_t exsig) const;

    ham_comp_t get_energy(const FrmBosHam& ham) const {
        return get_energy(ham, exsig::ex_1101) + get_energy(ham, exsig::ex_1110);
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
        return get_energy(ham.m_frm) + get_energy(ham.m_frmbos) + get_energy(ham.m_bos);
    }

private:
    void load_fn(const hdf5::NodeReader& /*parent*/) override {

    }

    void save_fn(const hdf5::NodeWriter& parent) override;
};

#endif //M7_RDMS_H
