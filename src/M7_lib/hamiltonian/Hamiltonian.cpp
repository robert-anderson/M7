//
// Created by Robert J. Anderson on 28/07/2021.
//

#include "Hamiltonian.h"
#include "M7_lib/hamiltonian/frm/GeneralFrmHam.h"
#include "M7_lib/hamiltonian/frmbos/GeneralLadderHam.h"
#include "M7_lib/hamiltonian/bos/InteractingBoseGasBosHam.h"

Hamiltonian::Hamiltonian(opt_pair_t opts): Hamiltonian(HamiltonianTerms(opts), nullptr, nullptr, nullptr){}


Hamiltonian::Hamiltonian(const FrmHam &frm): Hamiltonian({}, &frm, nullptr, nullptr){}


bool Hamiltonian::complex_valued() const {
    return m_frm.m_complex_valued;
}
