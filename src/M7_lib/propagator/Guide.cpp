//
// Created by rja on 13/08/22.
//

#include "Guide.h"

ham_t guide::Wavefunction::overlap(const FrmOnv& mbf) const {
    return frm_overlap(mbf);
}

ham_t guide::Wavefunction::overlap(const BosOnv& mbf) const {
    return bos_overlap(mbf);
}

ham_t guide::Wavefunction::overlap(const FrmBosOnv& mbf) const {
    return frmbos_overlap(mbf);
}

guide::EnergyDependent::EnergyDependent(const Hamiltonian& h) : m_h(h){}

guide::GutzwillerLike::GutzwillerLike(const Hamiltonian& h, double fac) : EnergyDependent(h), m_fac(fac){}
