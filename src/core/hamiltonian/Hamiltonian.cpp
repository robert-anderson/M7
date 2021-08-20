//
// Created by rja on 28/07/2021.
//

#include "Hamiltonian.h"

Hamiltonian::Hamiltonian(std::string fname, bool spin_major, size_t nboson_max, defs::ham_comp_t boson_frequency,
                         defs::ham_comp_t boson_coupling) :
        m_nboson_max(nboson_max), m_frm(fname, spin_major),
        m_frmbos(m_frm.m_nsite, m_nboson_max, fname),
        m_bos(m_frm.m_nsite, m_nboson_max, fname) {
}

Hamiltonian::Hamiltonian(const fciqmc_config::Hamiltonian &opts) :
        Hamiltonian(opts.m_fcidump.m_path, opts.m_fcidump.m_spin_major,
                    opts.m_nboson_max, opts.m_boson_frequency, opts.m_boson_coupling) {}

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
