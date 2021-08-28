//
// Created by rja on 05/11/2020.
//

#ifndef M7_HAMILTONIAN_H
#define M7_HAMILTONIAN_H

#include <type_traits>
#include <src/defs.h>
#include "FermionHamiltonian.h"
#include "BosonHamiltonian.h"
#include "LadderHamiltonian.h"
#include "src/core/nd/NdArray.h"
#include "HamiltonianParts.h"
#include "SymmetryHelpers.h"

using namespace field;

/**
 * generalized Hamiltonian class for fermionic, bosonic, and fermion-boson coupled interactions
 */
struct Hamiltonian {
    /**
     * the maximum number of bosons permitted to occupy any mode
     */
    const size_t m_nboson_max;
    /**
     * purely fermionic number-conserving terms in the Hamiltonian for traditional electronic structure calculations
     */
    FermionHamiltonian m_frm;
    /**
     * hamiltonian encapsulating all terms involving a single boson creation or annihilation operator
     * i.e. ranksigs 0010, 0001, 1110, 1101
     */
    LadderHamiltonian m_ladder;
    /**
     * purely bosonic number-conserving terms in the Hamiltonian
     */
    BosonHamiltonian m_bos;

    Hamiltonian(std::string fname, std::string fname_eb, std::string fname_bos, bool spin_major, size_t nboson_max = 0);

    Hamiltonian(std::string fname, bool spin_major): Hamiltonian(fname, "", "", spin_major, 0ul){}

    explicit Hamiltonian(const fciqmc_config::Hamiltonian &opts);

    size_t nci() const;

    const size_t &nsite() const;

    const size_t &nelec() const;

    bool complex_valued() const;

    /*
     * pure fermion matrix elements
     */

    defs::ham_t get_element(const FrmOnv &onv, const conn::FrmOnv &conn) const {
        return m_frm.get_element(onv, conn);
    }

    defs::ham_t get_element(const FrmOnv &onv) const {
        return m_frm.get_element_0000(onv);
    }

    defs::ham_comp_t get_energy(const FrmOnv &onv) const {
        return m_frm.get_energy(onv);
    }

    /*
     * pure boson matrix elements
     */

    defs::ham_t get_element(const BosOnv &onv, const conn::BosOnv &conn) const {
        return m_bos.get_element(onv, conn);
    }

    defs::ham_t get_element(const BosOnv &onv) const {
        return 0.0;
    }

    defs::ham_comp_t get_energy(const BosOnv &onv) const {
        return 0.0;
    }

    /*
     * fermion-boson coupled matrix elements
     */

    defs::ham_t get_element(const FrmBosOnv &onv, const conn::FrmBosOnv &conn) const {
        defs::ham_t helement_frm = 0.0;
        defs::ham_t helement_bos = 0.0;
        if (!conn.m_bos.size()) helement_frm = m_frm.get_element(onv.m_frm, conn.m_frm);
        if (!conn.m_frm.size()) helement_bos = m_bos.get_element(onv.m_bos, conn.m_bos);
        defs::ham_t helement_frmbos = m_ladder.get_element(onv, conn);
        return helement_frm + helement_bos + helement_frmbos;
    }

    defs::ham_t get_element(const FrmBosOnv &onv) const {
        return m_frm.get_element(onv.m_frm)+m_bos.get_element(onv.m_bos);
    }

    defs::ham_comp_t get_energy(const FrmBosOnv &onv) const {
        return consts::real(get_element(onv));
    }

    void set_hf_mbf(FrmOnv &onv, int spin) const;

    void set_hf_mbf(FrmBosOnv &onv, int spin) const;

    void set_afm_mbf(FrmOnv &onv, bool alpha_first) const;

    void set_afm_mbf(FrmBosOnv &onv, bool alpha_first) const;
};

#endif //M7_HAMILTONIAN_H
