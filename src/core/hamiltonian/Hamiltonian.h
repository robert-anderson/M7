//
// Created by rja on 05/11/2020.
//

#ifndef M7_HAMILTONIAN_H
#define M7_HAMILTONIAN_H

#include <type_traits>
#include <src/defs.h>
#include "FrmHam.h"
#include "BosHam.h"
#include "LadderHam.h"
#include "src/core/nd/NdArray.h"
#include "HubbardFrmHam.h"
#include "GeneralFrmHam.h"
#include "HolsteinLadderHam.h"
#include "HolsteinBosHam.h"

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
    const BasisDims m_bd;

private:
    BasisDims make_bd() const;

    std::unique_ptr<FrmHam> make_frm(const fciqmc_config::FermionHamiltonian &opts) {
        if (opts.m_hubbard)
            return std::unique_ptr<FrmHam>(new HubbardFrmHam(opts));
        else if (defs::enable_fermions)
            return std::unique_ptr<FrmHam>(new GeneralFrmHam(opts));
        return nullptr;
    }

    std::unique_ptr<LadderHam> make_ladder(const fciqmc_config::LadderHamiltonian &opts, size_t nsite) {
        if (opts.m_holstein_coupling) {
            auto nboson_max = opts.m_nboson_max.get();
            auto g = opts.m_holstein_coupling.get();
            return std::unique_ptr<LadderHam>(new HolsteinLadderHam(nsite, nboson_max, g));
        }
        return nullptr;
    }

    std::unique_ptr<BosHam> make_bos(const fciqmc_config::BosonHamiltonian &opts, size_t nsite) {
        if (opts.m_holstein_omega) {
            auto omega = opts.m_holstein_omega.get();
            return std::unique_ptr<BosHam>(new HolsteinBosHam(nsite, omega));
        }
        return std::unique_ptr<BosHam>(new BosHam(0ul, 0ul));
    }

    /*
    std::unique_ptr<LadderHamiltonian> make_ladder(const fciqmc_config::LadderHamiltonian &opts) {
        if (opts.m_holstein_coupling!=0.0)
            return std::unique_ptr<LadderHamiltonian>(new HubbardFrmHam(opts));
        else if (opts.m_ebdump)
            return std::unique_ptr<LadderHamiltonian>(new LadderHamiltonian(opts));
        else
            return nullptr;
    }
     */

public:

#if 0
    Hamiltonian(std::unique_ptr<FrmHam> &&frm,
                std::unique_ptr<LadderHam> &&ladder,
                std::unique_ptr<BosHam> &&bos) :
            m_frm(std::move(frm)), m_ladder(std::move(ladder)), m_bos(std::move(bos)),
            m_nboson_max(ladder ? ladder->m_nboson_max : 0ul), m_bd(make_bd()) {

    }

    Hamiltonian(std::string fname, std::string fname_eb, std::string fname_bos,
                bool spin_major, size_t nboson_max);

    Hamiltonian(std::string fname, bool spin_major) :
            Hamiltonian(fname, "", "", spin_major, 0ul) {}
#endif

    explicit Hamiltonian(const fciqmc_config::Hamiltonian &opts):
        m_frm(make_frm(opts.m_fermion)),
        m_ladder(make_ladder(opts.m_ladder, m_frm->m_nsite)),
        m_bos(make_bos(opts.m_boson, m_frm->m_nsite)),
        m_nboson_max(m_ladder ? m_ladder->m_nboson_max : 0ul), m_bd(make_bd()){
        if (!m_frm) log::info("Fermion Hamiltonian is disabled");
        if (defs::enable_bosons) {
            if (!m_ladder) log::info("Fermion-boson ladder Hamiltonian is disabled");
            if (!m_bos) log::info("Number-conserving boson Hamiltonian is disabled");
        }
    }

    size_t nci() const;

    size_t nelec() const;

    size_t nboson() const;

    /*
     * pure fermion coefficients
     */
    defs::ham_t get_coeff_1100(const size_t &i, const size_t &j) const {
        return m_frm->get_coeff_1100(i, j);
    }

    defs::ham_t get_coeff_2200(const size_t &i, const size_t &j, const size_t &k, const size_t &l) const {
        return m_frm->get_coeff_2200(i, j, k, l);
    }

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

    bool complex_valued() const {
        return m_frm->m_complex_valued;
    }
};

#endif //M7_HAMILTONIAN_H
