//
// Created by rja on 28/07/2021.
//

#include "Hamiltonian.h"
#include "GeneralFrmHam.h"

BasisDims Hamiltonian::make_bd() const {
    size_t nsite = 0ul;
    size_t nmode = 0ul;
    if (m_frm) nsite = m_frm->m_nsite;
    else if (m_ladder) return m_ladder->m_bd;
    if (m_bos) nmode = m_bos->m_nmode;
    return {nsite, nmode};
}

#if 0
Hamiltonian::Hamiltonian(std::string fname, std::string fname_eb, std::string fname_bos,
                         bool spin_major, size_t nboson_max):
        Hamiltonian(
                std::unique_ptr<FrmHam>(new GeneralFrmHam(fname, spin_major)),
                std::unique_ptr<LadderHam>(new LadderHam(fname_eb, nboson_max)),
                std::unique_ptr<BosHam>(new BosHam(fname_bos))){
    if (nboson_max && !defs::enable_bosons)
        log::warn("non-zero boson occupation specified with bosons disabled at compile time");
    if (nboson_max) {
        REQUIRE_EQ(m_ladder->m_bd.m_nsite, m_frm->m_nsite, "EBDUMP incompatible with FCIDUMP");
        REQUIRE_EQ(m_ladder->m_bd.m_nmode, m_bos->m_nmode, "EBDUMP incompatible with BOSDUMP");
    }
}
#endif

size_t Hamiltonian::nci() const {
    return m_frm->nci() * m_bos->nci();
}

size_t Hamiltonian::nelec() const {
    if (!m_frm) return 0ul;
    return m_frm->m_nelec;
}

size_t Hamiltonian::nboson() const {
    if (!m_bos) return 0ul;
    return m_bos->m_nboson;
}