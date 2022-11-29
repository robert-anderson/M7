//
// Created by Robert J. Anderson on 28/07/2021.
//

#include "Hamiltonian.h"
#include "M7_lib/hamiltonian/frm/GeneralFrmHam.h"
#include "M7_lib/hamiltonian/frmbos/GeneralLadderHam.h"
#include "M7_lib/hamiltonian/bos/InteractingBoseGasBosHam.h"

Hamiltonian::Hamiltonian(opt_pair_t opts): Hamiltonian(HamiltonianTerms(opts), nullptr, nullptr, nullptr){}

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