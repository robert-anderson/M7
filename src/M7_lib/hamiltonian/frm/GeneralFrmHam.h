//
// Created by anderson on 12/8/21.
//

#ifndef M7_GENERALFRMHAM_H
#define M7_GENERALFRMHAM_H

#include "FrmHam.h"

struct GeneralFrmHam : FrmHam {

    typedef Integrals_1e<defs::ham_t, defs::isym_1e> ints1_t;
    typedef Integrals_2e<defs::ham_t, defs::isym_2e> ints2_t;
    ints1_t m_int_1;
    ints2_t m_int_2;

    GeneralFrmHam(size_t nelec, size_t nsite, bool spin_resolved, int ms2_restrict, defs::inds site_irreps = {});

    GeneralFrmHam(const FcidumpHeader& header, bool spin_major, int ms2_restrict, int charge = 0);

    GeneralFrmHam(std::string fname, bool spin_major, int charge = 0):
            GeneralFrmHam(FcidumpHeader(fname), spin_major, charge){}

    explicit GeneralFrmHam(const fciqmc_config::FermionHamiltonian &opts) :
            GeneralFrmHam(opts.m_fcidump.m_path, opts.m_fcidump.m_spin_major, opts.m_charge) {}

    defs::ham_t get_coeff_1100(size_t i, size_t j) const override;

    defs::ham_t get_coeff_2200(size_t i, size_t j, size_t k, size_t l) const override;


    defs::ham_t get_element_0000(const field::FrmOnv &onv) const override;

    defs::ham_t get_element_1100(const field::FrmOnv &onv, const conn::FrmOnv &conn) const override;

    defs::ham_t get_element_2200(const field::FrmOnv &onv, const conn::FrmOnv &conn) const override;

    buffered::FrmOnv guess_reference(const int &spin_level) const;

    /*
    void add_excitgens() override {
        HamOpTerm::add_excitgens();
    }

    conn_foreach_ptr_list_t conn_foreach() override {
        conn_foreach_ptr_list_t list;
        if (m_kramers_attrs.m_conserving_singles) {
            //list.push_front(new conn_foreach::frm::General<>{})
        }
        return list;
    }
     */

private:
    /**
     * if electrons are disabled, the m_nsite members of the integrals will be set to zero so as
     * not to allocate unnecessary memory
     * @return
     * number of distinctly specified (spin-)orbitals in the Hamiltonian defintion
     */
    size_t norb_distinct() const {
        return (1ul + m_int_1.m_spin_res) * m_nsite;
    }
};


#endif //M7_GENERALFRMHAM_H
