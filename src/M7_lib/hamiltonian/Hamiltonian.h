//
// Created by Robert J. Anderson on 05/11/2020.
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
#include "M7_lib/hamiltonian/bos/NumOpBosHam.h"
#include "M7_lib/hamiltonian/bos/GeneralBosHam.h"
#include "M7_lib/hamiltonian/bos/HubbardBosHam.h"
#include "M7_lib/hamiltonian/bos/TcBosHam.h"
#include "M7_lib/hamiltonian/frm/HeisenbergFrmHam.h"
#include "M7_lib/hamiltonian/frm/SumFrmHam.h"
#include "M7_lib/hamiltonian/frm/SpinSquareFrmHam.h"
#include "M7_lib/hamiltonian/frm/TcFrmHam.h"
#include "M7_lib/hamiltonian/bos/InteractingBoseGasBosHam.h"
#include "M7_lib/hamiltonian/frmbos/GeneralLadderHam.h"
#include "M7_lib/util/Pointer.h"
#include "M7_lib/field/Mbf.h"


using namespace field;


/**
 * to employ polymorphism in the fermion, boson, and fermion-boson product terms of the hamiltonian, these must be
 * dynamically allocated and stored as pointers to the appropriate base classes. these are managed as unique_ptrs
 *
 * Precedence is always given to the FrmHam component, which never makes reference to the Bos or FrmBos components in
 * its creation
 *
 * FrmHam is set up first, so BosHam initialization can make read-only reference to it.
 * FrmBosHam is set up last, so its initialization can make read-only reference to both the FrmHam and BosHam.
 */
struct HamiltonianTerms {
    typedef HamOpTerm::InitOpts<conf::Hamiltonian> init_opts_t;
    /**
     * purely fermionic number-conserving terms in the Hamiltonian for traditional electronic structure calculations
     */
    std::shared_ptr<FrmHam> m_frm = nullptr;
    /**
     * purely bosonic number-conserving and non-conserving terms in the Hamiltonian
     */
    std::shared_ptr<BosHam> m_bos = nullptr;
    /**
     * hamiltonian encapsulating all terms involving products of fermion and boson operators
     */
    std::shared_ptr<FrmBosHam> m_frmbos = nullptr;
    /**
     * make the type of fermion Hamiltonian called for by the configuration, Then either return it directly, or combine
     * it with any modification specified in the options
     * @tparam ham_t
     *  FrmHam-derived class defining the Hamiltonian being created
     * @param opts
     *  configuration document section pertaining to the fermion hamiltonian
     * @return
     *  unique pointer to the polymorphic base class, FrmHam
     */
    template<typename ham_t>
    static std::shared_ptr<FrmHam> make_frm_modified(FrmHam::init_opts_t opts) {
        using namespace ptr::smart;
        static_assert(std::is_base_of<FrmHam, ham_t>::value, "template arg must be derived from FrmHam");
        ham_t bare_ham(opts);
        const auto spin_penalty_j = opts.m_ham.m_spin_penalty_j.m_value;
        if (spin_penalty_j != 0.0) {
            logging::info("Adding spin penalty to fermion Hamiltonian with J={}", spin_penalty_j);
            typedef SumFrmHam<ham_t, SpinSquareFrmHam> mod_ham_t;
            const FrmHam& base = bare_ham;
            const sys::frm::Sector sector(base.m_basis, base.electrons(opts.m_particles));
            return make_poly_shared<FrmHam, mod_ham_t>(std::move(bare_ham), SpinSquareFrmHam(sector), spin_penalty_j);
        }
        /*
         * if the configuration doc calls for no modification, just return the bare hamiltonian
         */
        return make_poly_shared<FrmHam, ham_t>(std::move(bare_ham));
    }

    std::shared_ptr<FrmHam> make_frm(FrmHam::init_opts_t opts) const;

    std::shared_ptr<BosHam> make_bos(BosHam::init_opts_t opts) const;

    std::shared_ptr<FrmBosHam> make_frmbos(FrmBosHam::init_opts_t opts) const;

    explicit HamiltonianTerms(init_opts_t opts);

    HamiltonianTerms(): m_frm(new NullFrmHam), m_bos(new NullBosHam), m_frmbos(new NullFrmBosHam){}

};


/**
 * generalized Hamiltonian class for fermionic, bosonic, and fermion-boson coupled interactions
 */
class Hamiltonian {
    typedef HamiltonianTerms::init_opts_t init_opts_t;

    const HamiltonianTerms m_terms;

public:
    /*
     * term ptrs are always dereferencable, so they are exposed as public const refs to the term base classes:
     */
    const FrmHam &m_frm;
    const BosHam &m_bos;
    const FrmBosHam &m_frmbos;

    const sys::Basis m_basis;
    /**
     * true if the Hamiltonian describes a close quantum system in the bosonic sector
     */
    const bool m_boson_number_conserve;

private:
    /**
     * workspace for computing connections
     */
    mutable suite::Conns m_work_conn;

    bool boson_number_conserve() const {
        if (m_frmbos.m_contribs_1101.any_nonzero()) return false;
        if (m_frmbos.m_contribs_1110.any_nonzero()) return false;
        return true;
    }

    /**
     * basic ctor kept private. this enables the flexibility to allow the FrmHam, BosHam, and FrmBosHam terms to be
     * owned either by the m_terms member or externally
     * @param terms
     *  moving reference to the terms object. can be trivially constructed in the case of external ownership
     * @param frm
     *  nullptr if m_terms.m_frm is to be dereferenced, else this points to an externally allocated FrmHam
     * @param bos
     *  nullptr if m_terms.m_bos is to be dereferenced, else this points to an externally allocated BosHam
     * @param frmbos
     *  nullptr if m_terms.m_frmbos is to be dereferenced, else this points to an externally allocated FrmBosHam
     */
    explicit Hamiltonian(HamiltonianTerms&& terms, const FrmHam* frm, const BosHam* bos, const FrmBosHam* frmbos);

    static void require_non_null(const HamOpTerm* ptr) {
        REQUIRE_TRUE(ptr, "pointer to externally-allocated term Hamiltonian must be non-null");
    }

public:

    /**
     * initialize based on the contents of a configuration document
     * @param opts
     *  pair of Sections from the configuration document: hamiltonian and basis
     */
    explicit Hamiltonian(init_opts_t opts);

    /*
     * ctors for initialization using externally-owned term objects (m_terms is initialized to nulls and referred to for
     * the terms
     */
    explicit Hamiltonian(const FrmHam* frm_ham);
    explicit Hamiltonian(const BosHam* bos_ham);
    explicit Hamiltonian(const FrmBosHam* bos_ham);
    Hamiltonian(const FrmHam* frm_ham, const FrmBosHam* frmbos_ham, const BosHam *bos_ham);


    /*
     * pure fermion matrix elements
     */

    ham_t get_element(const FrmOnv &onv, const conn::FrmOnv &conn) const {
        return m_frm.get_element(onv, conn);
    }

    ham_t get_element(const FrmOnv &onv) const {
        return m_frm.get_element_0000(onv);
    }

    ham_comp_t get_energy(const FrmOnv &onv) const {
        return m_frm.get_energy(onv);
    }

    /*
     * pure boson matrix elements
     */

    ham_t get_element(const BosOnv &onv, const conn::BosOnv &conn) const {
        return m_bos.get_element(onv, conn);
    }

    ham_t get_element(const BosOnv &onv) const {
        return m_bos.get_element(onv);
    }

    ham_comp_t get_energy(const BosOnv &onv) const {
        return m_bos.get_energy(onv);
    }

    /*
     * fermion-boson coupled matrix elements
     */

    ham_t get_element(const FrmBosOnv &onv, const conn::FrmBosOnv &conn) const {
        ham_t helement_frm = 0.0;
        ham_t helement_bos = 0.0;
        if (!conn.m_bos.size()) helement_frm = m_frm.get_element(onv.m_frm, conn.m_frm);
        if (!conn.m_frm.size()) helement_bos = m_bos.get_element(onv.m_bos, conn.m_bos);
        ham_t helement_ladder = m_frmbos.get_element(onv, conn);
        return helement_frm + helement_bos + helement_ladder;
    }

    ham_t get_element(const FrmBosOnv &onv) const {
        return get_element(onv.m_frm) + get_element(onv.m_bos);
    }

    ham_comp_t get_energy(const FrmBosOnv &onv) const {
        return arith::real(get_element(onv));
    }

    /*
     * convenience methods for matrix elements directly from bra and ket, using working connection object
     */
    template<typename mbf_t>
    ham_t get_element(const mbf_t &src, const mbf_t &dst) const {
        auto &conn = m_work_conn[src];
        // don't go to the expense of forming the connection if the excitation signature would be out of range
        auto exsig = mbf::exsig(src, dst);
        if (exsig == opsig::c_invalid) return 0.0;
        if (exsig.conserves_nfrm()) return 0.0;
        if (exsig.nfrm_cre() > 3ul) return 0.0;
        conn.connect(src, dst);
        return get_element(src, conn);
    }


    template<typename mbf_t>
    ham_t connected(const mbf_t &src, const mbf_t &dst) const {
        return ham::is_significant(get_element(src, dst));
    }

    bool complex_valued() const;

    // TODO: rename
    sys::Particles default_particles(uint_t nelec=0ul, int ms2=sys::frm::c_undefined_ms2, uint_t nboson=0ul) const;

    sys::Particles default_particles(const conf::Particles &opts) const;

    bool is_hermitian() const;

    /**
     * determine whether this Hamiltonian has a Brillouin theorem with respect to the given MBF
     */
    template<typename mbf_t>
    bool has_brillouin_theorem(const mbf_t&) const {
        return false;
    }

    /**
     * @param onv
     *  candidate Hartree-Fock determinant
     * @return
     *  true if there are single excitations from ONV with significant H matrix element
     */
    bool has_brillouin_theorem(const field::FrmOnv& onv) const;
};

#endif //M7_HAMILTONIAN_H
