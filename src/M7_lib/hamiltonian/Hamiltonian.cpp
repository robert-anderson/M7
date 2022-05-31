//
// Created by Robert J. Anderson on 28/07/2021.
//

#include "Hamiltonian.h"
#include "M7_lib/hamiltonian/frm/GeneralFrmHam.h"
#include "M7_lib/hamiltonian/frmbos/GeneralLadderHam.h"
#include "M7_lib/hamiltonian/bos/InteractingBoseGasBosHam.h"

Hamiltonian::Hamiltonian(opt_pair_t opts): Hamiltonian(HamiltonianTerms(opts), nullptr, nullptr, nullptr){}

Hamiltonian::Hamiltonian(const FrmHam *frm): Hamiltonian({}, frm, nullptr, nullptr){
    require_non_null(frm);
}

Hamiltonian::Hamiltonian(const BosHam *bos) : Hamiltonian({}, nullptr, bos, nullptr){
    require_non_null(bos);
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
