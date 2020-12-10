/**
 * @file
 * @author  Robert J. Anderson <robert.anderson@kcl.ac.uk>
 *
 * @section DESCRIPTION
 *
 * The fermion-boson Hamiltonian is accessed via instances of the (antisymmetric) fermion-boson
 * ONV connection classes.
 *
 * get_element_0: sum of fermionic and bosonic diagonal H matrix elements
 * get_element_1: purely fermionic 1-body H matrix element
 * get_element_2: purely fermionic 2-body H matrix element
 * get_element_01: H matrix element of the coupling between boson modes and (diagonal)
 *      fermionic density
 * // get_element_11: H matrix element of the coupling between boson modes and one-body fermionic
 *      excitations (not implemented for initial application: Hubbard--Holstein model).
 *
 *
 * get_energy: as in the purely fermionic case, returns the real part of get_element_0
 * get_element: delegating method where excitation levels must be queried prior to dispatch of one
 *      of the above methods.
 *
 */

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


    FermiBosHamiltonian(const Options &opts):
    FermiBosHamiltonian(opts.fcidump_path, opts.fcidump_spin_major, opts.nboson_max,
                        opts.boson_coupling, opts.boson_frequency){}

    const BosonCouplings &bc() const {
        return m_boson_couplings;
    }

    defs::ham_t get_element_0(const conn::Antisym<1> &afbconn) const {
        ASSERT(!afbconn);
        return FermionHamiltonian::get_element_0(afbconn) + m_boson_couplings.get_element_0(afbconn);
    }

    defs::ham_t get_element_0(const views::FbOnv &onv) const {
        return FermionHamiltonian::get_element_0(onv.m_fonv) + m_boson_couplings.get_element_0(onv.m_bonv);
    }

    using FermionHamiltonian::get_element_1;
    defs::ham_t get_element_1(const conn::Antisym<1> &afbconn) const {
        ASSERT(!afbconn.m_bonvconn);
        return FermionHamiltonian::get_element_1(afbconn);
    }

    using FermionHamiltonian::get_element_2;
    defs::ham_t get_element_2(const conn::Antisym<1> &afbconn) const {
        ASSERT(!afbconn.m_bonvconn);
        return FermionHamiltonian::get_element_2(afbconn);
    }

    defs::ham_t get_element_01(const conn::Antisym<1> &afbconn) const {
        ASSERT(!afbconn.nexcit());
        return m_boson_couplings.get_element_1(afbconn);
    }

    defs::ham_comp_t get_energy(const views::FbOnv &onv) const {
        return consts::real(get_element_0(onv));
    }

    defs::ham_t get_element(const conn::Antisym<1> &afbconn) const {
        if (!afbconn) return get_element_0(afbconn);
        else if (!afbconn.m_bonvconn) {
            if (afbconn.nexcit() == 1)
                return FermionHamiltonian::get_element_1(afbconn);
            else if (afbconn.nexcit() == 2)
                return FermionHamiltonian::get_element_2(afbconn);
        }
        else if (!afbconn.nexcit()) {
            if (afbconn.m_bonvconn.nexcit()==1)
                return get_element_01(afbconn);
        }
        return 0.0;
    }

    defs::ham_t get_element(const views::FbOnv &bra, const views::FbOnv &ket) const {
        return get_element(conn::Antisym<1>(ket, bra));
    }

    elements::FermiBosOnv guess_reference(const int &spin_level) const {
        elements::FermiBosOnv tmp(m_nsite);
        tmp.m_fonv = FermionHamiltonian::guess_reference(spin_level);
        return tmp;
    }


//
//    defs::ham_t get_element_0(const conn::AsFermiBosOnv &afbconn) const {
//        return FermionHamiltonian::get_element_1(afbconn.m_aconn);
//    }
//
//    defs::ham_t get_element_1(const conn::AsFermiBosOnv &afbconn) const {
//        return FermionHamiltonian::get_element_1(afbconn.m_aconn);
//    }
//
//    defs::ham_t get_element_2(const conn::AsFermiBosOnv &afbconn) const {
//        return FermionHamiltonian::get_element_2(afbconn.m_aconn);
//    }
//
//    defs::ham_t get_element_2(const conn::AsFermionOnv &aconn) const {
//        return FermionHamiltonian::get_element_2(aconn);
//    }
//
//    defs::ham_t get_element_2(const size_t &i, const size_t &j, const size_t &k, const size_t &l) const {
//        return FermionHamiltonian::get_element_2(i, j, k, l);
//    }
//
//    defs::ham_t get_element(const conn::AsFermiBosOnv &afbconn) const {
//        defs::ham_t res = m_boson_couplings.get_element(afbconn.m_aconn, afbconn.m_bonvconn);
//        if (afbconn.m_bonvconn.nchanged_mode() == 0) res += FermionHamiltonian::get_element(afbconn.m_aconn);
//        return res;
//    }
//
//    defs::ham_t get_element(const views::FermiBosOnv &bra, const views::FermiBosOnv &ket) const {
//        return get_element(conn::AsFermiBosOnv(ket, bra));
//    }

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
