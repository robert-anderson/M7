//
// Created by rja on 28/07/2021.
//

#include "Hamiltonian.h"
#include "GeneralFrmHam.h"

BasisDims Hamiltonian::make_bd() const {
    if (m_ladder->enabled()) return m_ladder->m_bd;
    return {m_frm->m_nsite, m_bos->m_nmode};
}

size_t Hamiltonian::nci() const {
    return m_frm->nci() * m_bos->nci();
}

size_t Hamiltonian::nelec() const {
    return m_frm->m_nelec;
}

size_t Hamiltonian::nboson() const {
    return m_bos->m_nboson;
}
