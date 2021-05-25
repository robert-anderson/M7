//
// Created by rja on 05/11/2020.
//

#include "FermiBosHamiltonian.h"

FermiBosHamiltonian::FermiBosHamiltonian(std::string fname, bool spin_major, size_t nboson_cutoff, defs::ham_t v,
                                         defs::ham_t omega) :
        FermionHamiltonian(fname, spin_major),
        m_boson_couplings(nsite(), nboson_cutoff, v, omega) {}

FermiBosHamiltonian::FermiBosHamiltonian(const Options &opts) :
        FermiBosHamiltonian(opts.fcidump_path, opts.fcidump_spin_major, opts.nboson_max,
                            opts.boson_coupling, opts.boson_frequency){}

const BosonCouplings &FermiBosHamiltonian::bc() const {
    return m_boson_couplings;
}

defs::ham_t FermiBosHamiltonian::get_element_0(const conn::Antisym<1> &afbconn) const {
    ASSERT(!afbconn);
    return FermionHamiltonian::get_element_0(afbconn) + m_boson_couplings.get_element_0(afbconn);
}

defs::ham_t FermiBosHamiltonian::get_element_1(const conn::Antisym<1> &afbconn) const {
    ASSERT(!afbconn.m_bonvconn);
    return FermionHamiltonian::get_element_1(afbconn);
}

defs::ham_t FermiBosHamiltonian::get_element_2(const conn::Antisym<1> &afbconn) const {
    ASSERT(!afbconn.m_bonvconn);
    return FermionHamiltonian::get_element_2(afbconn);
}

defs::ham_t FermiBosHamiltonian::get_element_01(const conn::Antisym<1> &afbconn) const {
    ASSERT(!afbconn.nexcit());
    return m_boson_couplings.get_element_1(afbconn);
}

defs::ham_t FermiBosHamiltonian::get_element(const fields::Onv<1> &bra, const fields::Onv<1> &ket) const {
    return get_element(conn::Antisym<1>(ket, bra));
}

defs::ham_t FermiBosHamiltonian::get_element_0(const fields::Onv<1> &onv) const {
    return FermionHamiltonian::get_element_0(onv.m_frm) + m_boson_couplings.get_element_0(onv.m_bos);
}
