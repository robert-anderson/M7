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
#include "M7_lib/hamiltonian/frm/HeisenbergFrmHam.h"
#include "M7_lib/hamiltonian/frm/SumFrmHam.h"
#include "M7_lib/hamiltonian/frm/SpinSquareFrmHam.h"
#include "M7_lib/hamiltonian/bos/InteractingBoseGasBosHam.h"
#include "M7_lib/hamiltonian/frmbos/GeneralLadderHam.h"


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
    typedef HamOpTerm::OptPair<conf::Hamiltonian> opt_pair_t;
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
    static std::unique_ptr<FrmHam> make_frm(FrmHam::opt_pair_t opts) {
        static_assert(std::is_base_of<FrmHam, ham_t>::value, "template arg must be derived from FrmHam");
        return std::unique_ptr<FrmHam>(new ham_t(opts));
        //auto j = opts.m_spin_penalty_j.get();
        /*
         * if the scalar of the spin square operator is zero, just return the bare hamiltonian
         */
        //if (j==0.0) return std::unique_ptr<FrmHam>(new ham_t(opts));
        //TODO: pass basis / hilbert config section
#if 0
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
#endif
    }

    std::unique_ptr<FrmHam> make_frm(FrmHam::opt_pair_t opts) {
        if (opts.m_ham.m_hubbard.enabled())
            return make_frm<HubbardFrmHam>(opts);
        else if (opts.m_ham.m_heisenberg.enabled())
            return make_frm<HeisenbergFrmHam>(opts);
        else if (opts.m_ham.m_fcidump.enabled())
            return make_frm<GeneralFrmHam>(opts);
        return std::unique_ptr<FrmHam>(new NullFrmHam);
    }

    std::unique_ptr<FrmBosHam> make_frmbos(FrmBosHam::opt_pair_t opts) {
        if (opts.m_ham.m_holstein_coupling) {
            const auto g = opts.m_ham.m_holstein_coupling.get();
            return std::unique_ptr<FrmBosHam>(new HolsteinLadderHam(m_frm->m_basis, g, opts.m_basis.m_bos_occ_cutoff));
        }
        else if (opts.m_ham.m_ebdump.enabled()) {
            return std::unique_ptr<FrmBosHam>(new GeneralLadderHam(opts));
        }
        return std::unique_ptr<FrmBosHam>(new NullFrmBosHam);
    }

    std::unique_ptr<BosHam> make_bos(BosHam::opt_pair_t opts){
        if (opts.m_ham.m_num_op_weight) {
            const size_t nsite = m_frm->m_basis.m_nsite;
            const sys::bos::Basis basis(nsite);
            const auto omega = opts.m_ham.m_num_op_weight.get();
            return std::unique_ptr<BosHam>(new NumOpBosHam(basis, omega));
        }
        else if (opts.m_ham.m_interacting_bose_gas.enabled())
            return std::unique_ptr<BosHam>(new InteractingBoseGasBosHam(opts));
        else if (opts.m_ham.m_bosdump.enabled())
            return std::unique_ptr<BosHam>(new GeneralBosHam(opts));
        return std::unique_ptr<BosHam>(new NullBosHam);
    }

    HamiltonianTerms(opt_pair_t opts):
        m_frm(make_frm({opts.m_ham.m_fermion, opts.m_basis})),
        m_bos(make_bos({opts.m_ham.m_boson, opts.m_basis})),
        m_frmbos(make_frmbos({opts.m_ham.m_ladder, opts.m_basis})){
        if (*m_frm && *m_frmbos)
            REQUIRE_TRUE(m_frm->m_basis==m_frmbos->m_basis.m_frm, "incompatible fermion basis definitions");
        if (*m_bos && *m_frmbos)
            REQUIRE_TRUE(m_bos->m_basis==m_frmbos->m_basis.m_bos, "incompatible boson basis definitions");
    }
};


/**
 * generalized Hamiltonian class for fermionic, bosonic, and fermion-boson coupled interactions
 */
class Hamiltonian {
    typedef HamiltonianTerms::opt_pair_t opt_pair_t;

    const HamiltonianTerms m_terms;

public:
    /*
     * term ptrs are always dereferencable, so they are exposed as public const refs to the term base classes:
     */
    const FrmHam &m_frm;
    const BosHam &m_bos;
    const FrmBosHam &m_frmbos;
    /**
     * true if the Hamiltonian describes a close quantum system in the bosonic sector
     */
    const bool m_boson_number_conserve;

private:

    mutable suite::Conns m_work_conn;

    bool boson_number_conserve() const {
        if (m_frmbos.m_contribs_1101.any_nonzero()) return true;
        if (m_frmbos.m_contribs_1110.any_nonzero()) return true;
        return false;
    }

public:

    explicit Hamiltonian(opt_pair_t opts);

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
        auto &conn = m_work_conn[src];
        conn.connect(src, dst);
        return get_element(src, conn);
    }

    bool complex_valued() const;

    /**
     * due the FCIDUMP format and default definitions of model Hamiltoniansm, the number of particles in the system and
     * their constraints are determined partly by this class, and partly by the user in the configuration document, this
     * function decides those parameters and raises errors when faced with any inconsistencies.
     * @param opts
     *  configuration document section
     * @return
     *  electron and boson number data
     */
//    sys::Particles get_quanta(const conf::Hamiltonian &opts) const {
//        // TODO: move config opts around
//        sys::frm::Electrons elecs(0ul);
//        sys::bos::Bosons bosons(0ul);
//        return {elecs, bosons};
//    }

    sys::Particles default_particles(const conf::Particles &opts) const {
        auto nelec = opts.m_nelec.get();
        // give precedence to nelec value given by FrmHam
        if (!nelec && m_frm) nelec = m_frm.default_nelec();
        if (!nelec && m_frmbos) nelec = m_frmbos.default_nelec();

        bool ms2_conserve = true;
        if (m_frm) ms2_conserve = m_frm.m_kramers_attrs.conserving();
        // currently only FrmHam can break Kramers symmetry (FrmBosHam always commutes with Sz)

        int ms2_value = sys::frm::Ms2::lowest_value(nelec);
        if (m_frm) ms2_value = m_frm.default_ms2()

        auto ms2 = (opts.m_ms2.get()==~0) ?
                   default_ms2(nelec) : sys::frm::Ms2(opts.m_ms2.get(), m_kramers_attrs.conserving());
        return {nelec, ms2};


        auto frm_particles = m_frm.default_particles(opts);
        auto frmbos_particles = m_frmbos.default_particles(opts);
        auto bos_particles = m_bos.default_particles(opts);


        auto nboson = opts.m_nboson.get();
        if (!nboson) nboson = m_bos.default_nboson();
        return {{nelec, ms2}, {nboson, m_boson_number_conserve}};
    }

    sys::Basis basis() const {
        return {m_frm.m_basis, m_bos.m_basis};
    }
    sys::Size basis_size() const {
        return {m_frm.m_basis.m_nsite, m_bos.m_basis.m_nmode};
    }
};

#endif //M7_HAMILTONIAN_H
