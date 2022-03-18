//
// Created by rja on 05/11/2020.
//

#ifndef M7_HAMILTONIAN_H
#define M7_HAMILTONIAN_H

#include <type_traits>
#include <defs.h>
#include <basis/Suites.h>
#include <nd/NdArray.h>

#include "FrmHam.h"
#include "BosHam.h"
#include "LadderHam.h"
#include "HubbardFrmHam.h"
#include "GeneralFrmHam.h"
#include "HolsteinLadderHam.h"
#include "HolsteinBosHam.h"
#include "GeneralBosHam.h"
#include "HeisenbergFrmHam.h"

using namespace field;

/**
 * generalized Hamiltonian class for fermionic, bosonic, and fermion-boson coupled interactions
 */
struct Hamiltonian {
    /**
     * purely fermionic number-conserving terms in the Hamiltonian for traditional electronic structure calculations
     */
    std::unique_ptr<FrmHam> m_frm;
    /**
     * hamiltonian encapsulating all terms involving a single boson creation or annihilation operator
     * i.e. ranksigs 0010, 0001, 1110, 1101
     */
    std::unique_ptr<LadderHam> m_ladder;
    /**
     * purely bosonic number-conserving terms in the Hamiltonian
     */
    std::unique_ptr<BosHam> m_bos;
    /**
     * the maximum number of bosons permitted to occupy any mode
     */
    const size_t m_nboson_max;
    /**
     * specifies number of fermion sites and boson modes defining the single-particle basis
     */
    const BasisData m_bd;

private:
    mutable suite::Conns m_work_conn;

    BasisData make_bd() const;

    std::unique_ptr<FrmHam> make_frm(const fciqmc_config::FermionHamiltonian &opts);

    std::unique_ptr<LadderHam> make_ladder(const fciqmc_config::LadderHamiltonian &opts, size_t nsite);

    std::unique_ptr<BosHam> make_bos(const fciqmc_config::BosonHamiltonian &opts, size_t nsite);

public:

    explicit Hamiltonian(const fciqmc_config::Hamiltonian &opts);

    size_t nci() const;

    size_t nelec() const;

    size_t nboson() const;

    /*
     * pure fermion matrix elements
     */

    defs::ham_t get_element(const FrmOnv &onv, const conn::FrmOnv &conn) const {
        return m_frm->get_element(onv, conn);
    }

    defs::ham_t get_element(const FrmOnv &onv) const {
        return m_frm->get_element_0000(onv);
    }

    defs::ham_comp_t get_energy(const FrmOnv &onv) const {
        return m_frm->get_energy(onv);
    }

    /*
     * pure boson matrix elements
     */

    defs::ham_t get_element(const BosOnv &onv, const conn::BosOnv &conn) const {
        return m_bos->get_element(onv, conn);
    }

    defs::ham_t get_element(const BosOnv &onv) const {
        return m_bos->get_element(onv);
    }

    defs::ham_comp_t get_energy(const BosOnv &onv) const {
        return m_bos->get_energy(onv);
    }

    /*
     * fermion-boson coupled matrix elements
     */

    defs::ham_t get_element(const FrmBosOnv &onv, const conn::FrmBosOnv &conn) const {
        defs::ham_t helement_frm = 0.0;
        defs::ham_t helement_bos = 0.0;
        if (!conn.m_bos.size()) helement_frm = m_frm->get_element(onv.m_frm, conn.m_frm);
        if (!conn.m_frm.size()) helement_bos = m_bos->get_element(onv.m_bos, conn.m_bos);
        defs::ham_t helement_ladder = m_ladder->get_element(onv, conn);
        return helement_frm + helement_bos + helement_ladder;
    }

    defs::ham_t get_element(const FrmBosOnv &onv) const {
        return get_element(onv.m_frm) + get_element(onv.m_bos);
    }

    defs::ham_comp_t get_energy(const FrmBosOnv &onv) const {
        return consts::real(get_element(onv));
    }

    /*
     * convenience methods for matrix elements directly from bra and ket, using working connection object
     */
    template<typename mbf_t>
    defs::ham_t get_element(const mbf_t &src, const mbf_t &dst) const {
        auto& conn = m_work_conn[src];
        conn.connect(src, dst);
        return get_element(src, conn);
    }

    bool complex_valued() const;
};

#endif //M7_HAMILTONIAN_H
