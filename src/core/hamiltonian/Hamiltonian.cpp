//
// Created by rja on 28/07/2021.
//

#include "Hamiltonian.h"

Hamiltonian::Hamiltonian(std::string fname, std::string fname_eb, std::string fname_bos, bool spin_major, size_t nboson_max):
        m_nboson_max(nboson_max), m_frm(fname, spin_major),
        m_frmbos(m_frm.m_nsite, m_nboson_max, fname_eb),
        m_bos(m_frm.m_nsite, m_nboson_max, fname_bos) {
    if (nboson_max && !defs::enable_bosons)
        log::warn("non-zero boson occupation specified with bosons disabled at compile time");
}

Hamiltonian::Hamiltonian(const fciqmc_config::Hamiltonian &opts) :
        Hamiltonian(opts.m_fcidump.m_path, opts.m_fcidump.m_eb_path, opts.m_fcidump.m_bos_path, opts.m_fcidump.m_spin_major,
                    opts.m_nboson_max) {}

size_t Hamiltonian::nci() const {
    return m_frm.nci() * m_bos.nci();
}

const size_t &Hamiltonian::nsite() const {
    return m_frm.m_nsite;
}

const size_t &Hamiltonian::nelec() const {
    return m_frm.m_nelec;
}

bool Hamiltonian::complex_valued() const {
    return m_frm.m_complex_valued;
}
