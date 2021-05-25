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

#include <src/core/table/BufferedFields.h>
#include "FermionHamiltonian.h"
#include "BosonCouplings.h"


struct FermiBosHamiltonian : public FermionHamiltonian {
    BosonCouplings m_boson_couplings;
    const std::unique_ptr<ham_sym_helpers::FermiBos> m_bc_sym_helper;
public:
    FermiBosHamiltonian(std::string fname, bool spin_major,
                        size_t nboson_cutoff, defs::ham_t v, defs::ham_t omega);


    FermiBosHamiltonian(const Options &opts);

    const BosonCouplings &bc() const;

    defs::ham_t get_element_0(const conn::Antisym<1> &afbconn) const;

    defs::ham_t get_element_0(const fields::Onv<1> &onv) const;

    using FermionHamiltonian::get_element_1;
    defs::ham_t get_element_1(const conn::Antisym<1> &afbconn) const;

    using FermionHamiltonian::get_element_2;
    defs::ham_t get_element_2(const conn::Antisym<1> &afbconn) const;

    defs::ham_t get_element_01(const conn::Antisym<1> &afbconn) const;

    defs::ham_t get_element(const fields::Onv<1> &bra, const fields::Onv<1> &ket) const;

private:
    defs::ham_t get_element_tag(const conn::Antisym<1> &aconn, dispatch_utils::BoolTag<false> sym_opts) const {
        if (!aconn) return get_element_0(aconn);
        else if (!aconn.m_bonvconn) {
            if (aconn.nexcit() == 1)
                return FermionHamiltonian::get_element_1(aconn);
            else if (aconn.nexcit() == 2)
                return FermionHamiltonian::get_element_2(aconn);
        }
        else if (!aconn.nexcit()) {
            if (aconn.m_bonvconn.nexcit()==1)
                return get_element_01(aconn);
        }
        return 0.0;
    }

    defs::ham_t get_element_tag(const conn::Antisym<1> &aconn, dispatch_utils::BoolTag<true> sym_opts) const {
        return m_bc_sym_helper->get_element(aconn);
    }

    defs::ham_comp_t get_energy_tag(const fields::Onv<1> &onv, dispatch_utils::BoolTag<false> sym_opts) const {
        return consts::real(get_element_0(onv));
    }
    defs::ham_comp_t get_energy_tag(const fields::Onv<1> &onv, dispatch_utils::BoolTag<true> sym_opts) const {
        return m_bc_sym_helper->get_energy(onv);
    }

public:

    defs::ham_t get_element(const conn::Antisym<1> &aconn) const {
        return get_element_tag(aconn, dispatch_utils::BoolTag<defs::enable_optim_for_lattice_ham>());
    }

    defs::ham_comp_t get_energy(const fields::Onv<1> &onv) const{
        return get_energy_tag(onv, dispatch_utils::BoolTag<defs::enable_optim_for_lattice_ham>());
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

    void foreach_connection(const fields::Onv<1> &src_onv, const ham_sym_helpers::FermiBos::body_fn_t &body,
                            bool get_h, bool h_nonzero_only, bool include_diagonal) const {
        m_bc_sym_helper->foreach_connection(src_onv, body, get_h, h_nonzero_only, include_diagonal);
    }

    void set_hf_onv(fields::Onv<1>& onv, int spin) const {
        FermionHamiltonian::set_hf_onv(onv.m_frm, spin);
    }
};


#endif //M7_BOSONCOUPLINGS_H
