//
// Created by Robert J. Anderson on 28/07/2021.
//

#include "Hamiltonian.h"
#include "M7_lib/hamiltonian/frm/GeneralFrmHam.h"
#include "M7_lib/hamiltonian/frmbos/GeneralLadderHam.h"
#include "M7_lib/hamiltonian/bos/InteractingBoseGasBosHam.h"

Hamiltonian::Hamiltonian(opt_pair_t opts):
        m_terms(opts), m_frm(*m_terms.m_frm), m_bos(*m_terms.m_bos), m_frmbos(*m_terms.m_frmbos),
        m_basis(sys::frm::Basis(m_frmbos.m_basis.m_frm, m_frm.m_basis),
                sys::bos::Basis(m_frmbos.m_basis.m_bos, m_bos.m_basis)), m_work_conn(m_basis.size()){
    REQUIRE_TRUE(m_basis, "No system defined");
    if (m_frm.disabled()) log::info("Fermion Hamiltonian is disabled");
    if (defs::enable_bosons) {
        if (m_frmbos.disabled()) log::info("Fermion-boson ladder Hamiltonian is disabled");
        if (m_bos.disabled()) log::info("Number-conserving boson Hamiltonian is disabled");
    }
}

bool Hamiltonian::complex_valued() const {
    return m_frm.m_complex_valued;
}