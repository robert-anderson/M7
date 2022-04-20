//
// Created by rja on 28/07/2021.
//

#include "Hamiltonian.h"
#include "M7_lib/hamiltonian/frm/GeneralFrmHam.h"
#include "M7_lib/hamiltonian/frmbos/GeneralLadderHam.h"
#include "M7_lib/hamiltonian/bos/InteractingBoseGasBosHam.h"

Hamiltonian::Hamiltonian(const fciqmc_config::Hamiltonian &opts) :
        m_terms(opts), m_frm(*m_terms.m_frm), m_bos(*m_terms.m_bos), m_frmbos(*m_terms.m_frmbos),
        m_bd(m_frmbos.enabled() ? m_frmbos.m_bd : BasisData(m_frm.m_bd, m_bos.m_bd)), m_work_conn(m_bd){
    REQUIRE_TRUE(m_bd.m_frm.m_nsite || m_bd.m_bos.m_nmode, "No system defined");
    if (m_frm.disabled()) log::info("Fermion Hamiltonian is disabled");
    if (defs::enable_bosons) {
        if (m_frmbos.disabled()) log::info("Fermion-boson ladder Hamiltonian is disabled");
        if (m_bos.disabled()) log::info("Number-conserving boson Hamiltonian is disabled");
    }
}

size_t Hamiltonian::nelec() const {
    return m_frm.m_hd.m_nelec;
}

size_t Hamiltonian::nboson() const {
    return m_bos.m_hd.m_nboson;
}

bool Hamiltonian::complex_valued() const {
    return m_frm.m_complex_valued;
}