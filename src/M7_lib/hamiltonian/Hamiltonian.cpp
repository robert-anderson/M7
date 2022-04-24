//
// Created by rja on 28/07/2021.
//

#include "Hamiltonian.h"
#include "M7_lib/hamiltonian/frm/GeneralFrmHam.h"
#include "M7_lib/hamiltonian/frmbos/GeneralLadderHam.h"
#include "M7_lib/hamiltonian/bos/InteractingBoseGasBosHam.h"

Hamiltonian::Hamiltonian(const fciqmc_config::Hamiltonian &opts) :
        m_terms(opts), m_frm(*m_terms.m_frm), m_bos(*m_terms.m_bos), m_frmbos(*m_terms.m_frmbos),
        m_hs(m_frmbos.m_hs, HilbertSpace{m_frm.m_hs, m_bos.m_hs}), m_work_conn(m_hs.m_extents){
    REQUIRE_TRUE(m_hs.m_frm || m_hs.m_bos, "No system defined");
    if (m_frm.disabled()) log::info("Fermion Hamiltonian is disabled");
    if (defs::enable_bosons) {
        if (m_frmbos.disabled()) log::info("Fermion-boson ladder Hamiltonian is disabled");
        if (m_bos.disabled()) log::info("Number-conserving boson Hamiltonian is disabled");
    }
}

bool Hamiltonian::complex_valued() const {
    return m_frm.m_complex_valued;
}