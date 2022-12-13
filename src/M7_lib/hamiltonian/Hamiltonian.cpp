//
// Created by Robert J. Anderson on 28/07/2021.
//

#include "Hamiltonian.h"
#include "M7_lib/hamiltonian/frm/GeneralFrmHam.h"
#include "M7_lib/hamiltonian/frmbos/GeneralLadderHam.h"
#include "M7_lib/hamiltonian/bos/InteractingBoseGasBosHam.h"

std::shared_ptr<FrmHam> HamiltonianTerms::make_frm(FrmHam::init_opts_t opts) const {
    using namespace ptr::smart;
    if (opts.m_ham.m_hubbard.m_enabled)
        return make_frm_modified<HubbardFrmHam>(opts);
    else if (opts.m_ham.m_heisenberg.m_enabled)
        return make_frm_modified<HeisenbergFrmHam>(opts);
    else if (opts.m_ham.m_fcidump.m_enabled)
        return make_frm_modified<GeneralFrmHam>(opts);
    return make_poly_shared<FrmHam, NullFrmHam>();
}

std::shared_ptr<BosHam> HamiltonianTerms::make_bos(BosHam::init_opts_t opts) const {
    using namespace ptr::smart;
    if (opts.m_ham.m_num_op_weight) {
        const uint_t nsite = m_frm->m_basis.m_nsite;
        const sys::bos::Basis basis(nsite, opts.m_basis.m_bos_occ_cutoff);
        const auto omega = opts.m_ham.m_num_op_weight.m_value;
        return make_poly_shared<BosHam, NumOpBosHam>(basis, omega);
    }
    else if (opts.m_ham.m_interacting_bose_gas.m_enabled)
        return make_poly_shared<BosHam, InteractingBoseGasBosHam>(opts);
    else if (opts.m_ham.m_hubbard.m_enabled)
        return make_poly_shared<BosHam, HubbardBosHam>(opts);
    else if (opts.m_ham.m_bosdump.m_enabled)
        return make_poly_shared<BosHam, GeneralBosHam>(opts);
    return make_poly_shared<BosHam, NullBosHam>();
}

std::shared_ptr<FrmBosHam> HamiltonianTerms::make_frmbos(FrmBosHam::init_opts_t opts) const {
    REQUIRE_TRUE(m_frm.get(), "fermion Hamiltonian unallocated");
    REQUIRE_TRUE(m_bos.get(), "boson Hamiltonian unallocated");
    const sys::Basis basis(m_frm->m_basis, m_bos->m_basis);

    using namespace ptr::smart;
    if (opts.m_ham.m_holstein_coupling.m_value != 0.0) {
        const auto g = opts.m_ham.m_holstein_coupling.m_value;
        return make_poly_shared<FrmBosHam, HolsteinLadderHam>(basis, g);
    }
    else if (opts.m_ham.m_ebdump.m_enabled) {
        return make_poly_shared<FrmBosHam, GeneralLadderHam>(basis, opts);
    }
    return make_poly_shared<FrmBosHam, NullFrmBosHam>();
}

HamiltonianTerms::HamiltonianTerms(HamiltonianTerms::init_opts_t opts) :
        m_frm(make_frm({opts.m_ham.m_fermion, opts.m_basis, opts.m_particles})),
        m_bos(make_bos({opts.m_ham.m_boson, opts.m_basis, opts.m_particles})),
        m_frmbos(make_frmbos({opts.m_ham.m_ladder, opts.m_basis, opts.m_particles})){
    if (*m_frm && *m_frmbos)
        REQUIRE_TRUE(m_frm->m_basis==m_frmbos->m_basis.m_frm, "incompatible fermion basis definitions");
    if (*m_bos && *m_frmbos)
        REQUIRE_TRUE(m_bos->m_basis==m_frmbos->m_basis.m_bos, "incompatible boson basis definitions");
}

Hamiltonian::Hamiltonian(HamiltonianTerms&& terms, const FrmHam* frm, const BosHam* bos, const FrmBosHam* frmbos) :
        m_terms(std::move(terms)), m_frm(frm ? *frm : *m_terms.m_frm), m_bos(bos ? *bos : *m_terms.m_bos),
        m_frmbos(frmbos ? *frmbos : *m_terms.m_frmbos),
        m_basis(frmbos ? m_frmbos.m_basis : sys::Basis(m_frm.m_basis, m_bos.m_basis)),
        m_boson_number_conserve(boson_number_conserve()), m_work_conn(m_basis.size()){
    REQUIRE_TRUE(m_basis, "No system defined");
    if (!m_frm) logging::info("Fermion Hamiltonian is disabled");
    if (c_enable_bosons) {
        if (!m_frmbos) logging::info("Fermion-boson ladder Hamiltonian is disabled");
        if (!m_bos) logging::info("Number-conserving boson Hamiltonian is disabled");
    }

    if (m_frm && m_frmbos) REQUIRE_EQ(m_frm.m_basis, m_frmbos.m_basis.m_frm,
                                      "Frm and FrmBos H terms do not have the same fermionic basis definition");
    if (m_bos && m_frmbos) REQUIRE_EQ(m_bos.m_basis, m_frmbos.m_basis.m_bos,
                                      "Bos and FrmBos H terms do not have the same bosonic basis definition");
    logging::info("Hamiltonian is {}hermitian", (is_hermitian() ? "" : "NON-"));
    if (m_frm.m_e_core != 0.0) logging::info(
        "Hamiltonian has a non-zero core energy {} which is included in all energy estimators", m_frm.m_e_core);
}

Hamiltonian::Hamiltonian(init_opts_t opts): Hamiltonian(HamiltonianTerms(opts), nullptr, nullptr, nullptr){}

Hamiltonian::Hamiltonian(const FrmHam *ham): Hamiltonian({}, ham, nullptr, nullptr){
    require_non_null(ham);
}

Hamiltonian::Hamiltonian(const BosHam *ham) : Hamiltonian({}, nullptr, ham, nullptr){
    require_non_null(ham);
}

Hamiltonian::Hamiltonian(const FrmBosHam *ham) : Hamiltonian({}, nullptr, nullptr, ham){
    require_non_null(ham);
}

Hamiltonian::Hamiltonian(const FrmHam *frm, const FrmBosHam *frmbos, const BosHam *bos) :
        Hamiltonian({}, frm, bos, frmbos){
    require_non_null(frm);
    require_non_null(frmbos);
    require_non_null(bos);
}


bool Hamiltonian::complex_valued() const {
    return m_frm.m_complex_valued;
}

sys::Particles Hamiltonian::default_particles(uint_t nelec, int ms2, uint_t nboson) const {
    // if nelec is not already set, give precedence to nelec value given by FrmHam
    if (!nelec && m_frm) nelec = m_frm.default_nelec();
    if (!nelec && m_frmbos) nelec = m_frmbos.default_nelec();

    bool ms2_conserve = true;
    if (m_frm) ms2_conserve = m_frm.m_kramers_attrs.conserving();
    // currently only FrmHam can break Kramers symmetry (FrmBosHam always commutes with Sz)

    if (m_frm && ms2==sys::frm::c_undefined_ms2) ms2 = m_frm.default_ms2_value();
    if (m_frmbos && ms2==sys::frm::c_undefined_ms2) ms2 = m_frmbos.default_ms2_value();
    if (ms2==sys::frm::c_undefined_ms2) {
        logging::info("2*Ms value not defined by configuration document or Hamiltonian, "
                      "defaulting to lowest valid positive value");
        ms2 = sys::frm::Ms2::lowest_value(nelec);
    }

    // give precedence to nboson value given by BosHam
    if (!nboson) nboson = m_bos.default_nboson();
    if (!nboson) nboson = m_frmbos.default_nboson();

    return {{nelec, {ms2, ms2_conserve}}, {nboson, m_boson_number_conserve}};
}

sys::Particles Hamiltonian::default_particles(const conf::Particles &opts) const {
    return default_particles(opts.m_nelec, opts.m_ms2, opts.m_nboson);
}

bool Hamiltonian::is_hermitian() const {
    return m_frm.is_hermitian() && m_frmbos.is_hermitian() && m_bos.is_hermitian();
}

bool Hamiltonian::has_brillouin_theorem(const FrmOnv &onv) const {
    conn::FrmOnv conn(onv);
    bool any_nonzero = false;
    auto fn = [&](){any_nonzero |= ham::is_significant(get_element(onv, conn));};
    if (m_frm.m_kramers_attrs.conserving()) {
        conn_foreach::frm::Ms2Conserve<1> foreach;
        foreach.loop_fn(conn, onv, fn);
    }
    else {
        conn_foreach::frm::General<1> foreach;
        foreach.loop_fn(conn, onv, fn);
    }
    return !any_nonzero;
}