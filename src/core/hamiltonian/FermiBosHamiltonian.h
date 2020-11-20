//
// Created by rja on 05/11/2020.
//

#ifndef M7_FERMIBOSHAMILTONIAN_H
#define M7_FERMIBOSHAMILTONIAN_H

#include "FermionHamiltonian.h"
#include "BosonCouplings.h"

class FermiBosHamiltonian : public FermionHamiltonian {
    BosonCouplings m_boson_couplings;
public:
    FermiBosHamiltonian(std::string fname, bool spin_major,
                        size_t nboson_cutoff, defs::ham_t v, defs::ham_t omega) :
            FermionHamiltonian(fname, spin_major),
            m_boson_couplings(nsite(), nboson_cutoff, v, omega) {}

    const BosonCouplings &bc() const {
        return m_boson_couplings;
    }

    defs::ham_t get_element(const conn::AsFermiBosOnv &afbconn) const {
        defs::ham_t res = m_boson_couplings.get_element(afbconn.m_aconn, afbconn.m_bonvconn);
        if (afbconn.m_bonvconn.nchanged_mode() == 0) res += FermionHamiltonian::get_element(afbconn.m_aconn);
        return res;
    }

    defs::ham_t get_element(const views::FermiBosOnv &bra, const views::FermiBosOnv &ket) const {
        return get_element(conn::AsFermiBosOnv(ket, bra));
    }

    size_t nci() const {
        return FermionHamiltonian::nci() * m_boson_couplings.nci();
    }

    const size_t &nmode() const {
        return m_boson_couplings.nmode();
    }

    const size_t &nboson_cutoff() const {
        return m_boson_couplings.nboson_cutoff();
    }
};


#endif //M7_BOSONCOUPLINGS_H
