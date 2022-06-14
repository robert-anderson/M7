//
// Created by Robert J. Anderson on 28/07/2021.
//

#include "Hamiltonian.h"
#include "M7_lib/hamiltonian/frm/GeneralFrmHam.h"
#include "M7_lib/hamiltonian/frmbos/GeneralLadderHam.h"
#include "M7_lib/hamiltonian/bos/InteractingBoseGasBosHam.h"

Hamiltonian::Hamiltonian(opt_pair_t opts): Hamiltonian(HamiltonianTerms(opts), nullptr, nullptr, nullptr){}

Hamiltonian::Hamiltonian(const FrmHam *ham): Hamiltonian({}, ham, nullptr, nullptr){
    require_non_null(ham);
}

Hamiltonian::Hamiltonian(const BosHam *ham) : Hamiltonian({}, nullptr, ham, nullptr){
    require_non_null(ham);
}

Hamiltonian::Hamiltonian(const FrmBosHam *ham) : Hamiltonian({}, nullptr, nullptr, ham){
    require_non_null(ham);
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
