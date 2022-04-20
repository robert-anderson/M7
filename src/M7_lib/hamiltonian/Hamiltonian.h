//
// Created by rja on 05/11/2020.
//

#ifndef M7_HAMILTONIAN_H
#define M7_HAMILTONIAN_H

#include <type_traits>
#include <M7_lib/defs.h>
#include <M7_lib/basis/Suites.h>
#include <M7_lib/nd/NdArray.h>

#include "M7_lib/hamiltonian/frm/FrmHam.h"
#include "M7_lib/hamiltonian/bos/BosHam.h"
#include "M7_lib/hamiltonian/frmbos/FrmBosHam.h"
#include "M7_lib/hamiltonian/frm/HubbardFrmHam.h"
#include "M7_lib/hamiltonian/frm/GeneralFrmHam.h"
#include "M7_lib/hamiltonian/frmbos/HolsteinLadderHam.h"
#include "M7_lib/hamiltonian/bos/HolsteinBosHam.h"
#include "M7_lib/hamiltonian/bos/GeneralBosHam.h"
#include "M7_lib/hamiltonian/frm/HeisenbergFrmHam.h"
#include "M7_lib/hamiltonian/frm/SumFrmHam.h"
#include "M7_lib/hamiltonian/frm/SpinSquareFrmHam.h"
#include "M7_lib/hamiltonian/bos/InteractingBoseGasBosHam.h"
#include "M7_lib/hamiltonian/frmbos/GeneralLadderHam.h"


using namespace field;


/**
 * to employ polymorphism in the fermion, boson, and fermion-boson product terms of the hamiltonian must be
 * dynamically allocated and stored as a pointer to the appropriate base class. these are managed as unique_ptrs
 *
 * FrmHam is set up first, so BosHam initialization can make read-only reference to it.
 * FrmBosHam is set up last, so its initialization can make read-only reference to both the FrmHam and BosHam.
 */
struct HamiltonianTerms {
    /**
     * purely fermionic number-conserving terms in the Hamiltonian for traditional electronic structure calculations
     */
    std::unique_ptr<FrmHam> m_frm = nullptr;
    /**
     * purely bosonic number-conserving and non-conserving terms in the Hamiltonian
     */
    std::unique_ptr<BosHam> m_bos = nullptr;
    /**
     * hamiltonian encapsulating all terms involving products of fermion and boson operators
     */
    std::unique_ptr<FrmBosHam> m_frmbos = nullptr;
    /**
     * make the type of fermion Hamiltonian called for by the configuration, Then either return it directly, or combine
     * it with the spin penalty if this modification is specified in the options
     * @tparam ham_t
     *  FrmHam-derived class defining the Hamiltonian being created
     * @param opts
     *  configuration document section pertaining to the fermion hamiltonian
     * @return
     *  unique pointer to the polymorphic base class, FrmHam
     */
    template<typename ham_t>
    std::unique_ptr<FrmHam> make_frm(const fciqmc_config::FermionHamiltonian &opts) {
        static_assert(std::is_base_of<FrmHam, ham_t>::value, "template arg must be derived from FrmHam");
        auto j = opts.m_spin_penalty_j.get();
        /*
         * if the scalar of the spin square operator is zero, just return the bare hamiltonian
         */
        if (j==0.0) return std::unique_ptr<FrmHam>(new ham_t(opts));
        /*
         * build the bare hamiltonian and let a spin square hamiltonian instance be created by read-only access to the
         * bare hamiltonian
         */
        ham_t bare_ham(opts);
        SpinSquareFrmHam spin_ham(bare_ham);
        /*
         * now the sum can be cheaply created by moving these two components
         */
        return std::unique_ptr<FrmHam>(new SumFrmHam<ham_t, SpinSquareFrmHam>(std::move(bare_ham), std::move(spin_ham), j));
    }

    std::unique_ptr<FrmHam> make_frm(const fciqmc_config::FermionHamiltonian &opts) {
        if (opts.m_hubbard.enabled())
            return make_frm<HubbardFrmHam>(opts);
        else if (opts.m_heisenberg.enabled())
            return make_frm<HeisenbergFrmHam>(opts);
        else if (opts.m_fcidump.enabled())
            return make_frm<GeneralFrmHam>(opts);
        return std::unique_ptr<FrmHam>(new NullFrmHam);
    }

    std::unique_ptr<BosHam> make_bos(const fciqmc_config::BosonHamiltonian &opts, const FrmHam* frm){
        if (opts.m_holstein_omega) {
            auto omega = opts.m_holstein_omega.get();
            return std::unique_ptr<BosHam>(new HolsteinBosHam(frm->m_bd.m_nsite, opts.m_nboson, omega));
        }
        else if (opts.m_interacting_bose_gas.enabled())
            return std::unique_ptr<BosHam>(new InteractingBoseGasBosHam(opts));
        else if (opts.m_bosdump.enabled())
            return std::unique_ptr<BosHam>(new GeneralBosHam(opts));
        return std::unique_ptr<BosHam>(new NullBosHam);
    }

    std::unique_ptr<FrmBosHam> make_frmbos(const fciqmc_config::FrmBosHamiltonian &opts,
                                           const FrmHam* frm, const BosHam* bos) {
        if (opts.m_holstein_coupling) {
            REQUIRE_TRUE(dynamic_cast<const HubbardFrmHam*>(frm),
                         "Holstein coupling requires Hubbard-type fermion Hamiltonian");
            auto nsite = frm->m_bd.m_nsite;
            auto g = opts.m_holstein_coupling.get();
            return std::unique_ptr<FrmBosHam>(new HolsteinLadderHam({nsite, nsite}, g));
        }
        else if (opts.m_ebdump.enabled()) return std::unique_ptr<FrmBosHam>(new GeneralLadderHam(opts));
        return std::unique_ptr<FrmBosHam>(new NullLadderHam);
    }

    HamiltonianTerms(const fciqmc_config::Hamiltonian &opts):
        m_frm(make_frm(opts.m_fermion)),
        m_bos(make_bos(opts.m_boson, m_frm.get())),
        m_frmbos(make_frmbos(opts.m_ladder, m_frm.get(), m_bos.get())){}
};


/**
 * generalized Hamiltonian class for fermionic, bosonic, and fermion-boson coupled interactions
 */
struct Hamiltonian {
    HamiltonianTerms m_terms;

public:
    /*
     * term ptrs are always dereferencable, so they are exposed as public const refs to the base classes:
     */
    const FrmHam& m_frm;
    const BosHam& m_bos;
    const FrmBosHam& m_frmbos;
    /**
     * specifies number of fermion sites and boson modes along with other attributes defining the single-particle basis
     */
    const BasisData m_bd;

private:
    mutable suite::Conns m_work_conn;

public:

    explicit Hamiltonian(const fciqmc_config::Hamiltonian &opts);

    size_t nelec() const;

    size_t nboson() const;

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
        return m_bos.get_element(onv);
    }

    defs::ham_comp_t get_energy(const BosOnv &onv) const {
        return m_bos.get_energy(onv);
    }

    /*
     * fermion-boson coupled matrix elements
     */

    defs::ham_t get_element(const FrmBosOnv &onv, const conn::FrmBosOnv &conn) const {
        defs::ham_t helement_frm = 0.0;
        defs::ham_t helement_bos = 0.0;
        if (!conn.m_bos.size()) helement_frm = m_frm.get_element(onv.m_frm, conn.m_frm);
        if (!conn.m_frm.size()) helement_bos = m_bos.get_element(onv.m_bos, conn.m_bos);
        defs::ham_t helement_ladder = m_frmbos.get_element(onv, conn);
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
