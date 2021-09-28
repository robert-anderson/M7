//
// Created by rja on 28/07/2021.
//

#include "Hamiltonian.h"

Hamiltonian::Hamiltonian(std::string fname, std::string fname_eb, std::string fname_bos, bool spin_major, size_t nboson_max):
        m_nboson_max(nboson_max), m_frm(fname, spin_major), m_bos(fname_bos, m_nboson_max),
        m_ladder(fname_eb, m_nboson_max), m_bd(m_frm.m_nsite, m_bos.m_nmode) {
    if (nboson_max && !defs::enable_bosons)
        log::warn("non-zero boson occupation specified with bosons disabled at compile time");
    if (nboson_max) {
        REQUIRE_EQ(m_ladder.m_bd.m_nsite, m_frm.m_nsite, "EBDUMP incompatible with FCIDUMP");
        REQUIRE_EQ(m_ladder.m_bd.m_nmode, m_bos.m_nmode, "EBDUMP incompatible with BOSDUMP");
    }
}

Hamiltonian::Hamiltonian(const fciqmc_config::Hamiltonian &opts) :
        Hamiltonian(opts.m_fcidump.m_path, opts.m_fcidump.m_eb_path, opts.m_fcidump.m_bos_path, opts.m_fcidump.m_spin_major,
                    opts.m_nboson_max) {}

size_t Hamiltonian::nci() const {
    return m_frm.nci() * m_bos.nci();
}

const size_t &Hamiltonian::nelec() const {
    return m_frm.m_nelec;
}

bool Hamiltonian::complex_valued() const {
    return m_frm.m_complex_valued;
}

void Hamiltonian::set_hf_mbf(FrmOnv &onv, int spin) const {
    m_frm.set_hf_mbf(onv, spin);
}

void Hamiltonian::set_hf_mbf(FrmBosOnv &onv, int spin) const {
    m_frm.set_hf_mbf(onv.m_frm, spin);
}

void Hamiltonian::set_afm_mbf(FrmOnv &onv, bool alpha_first) const {
    m_frm.set_afm_mbf(onv, alpha_first);
}

void Hamiltonian::set_afm_mbf(FrmBosOnv &onv, bool alpha_first) const {
    m_frm.set_afm_mbf(onv.m_frm, alpha_first);
}
