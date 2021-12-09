//
// Created by anderson on 12/8/21.
//

#ifndef M7_ABINITIOHAMILTONIAN_H
#define M7_ABINITIOHAMILTONIAN_H

#include "FermionHamiltonian.h"

struct AbinitioHamiltonian : FermionHamiltonian {

    typedef Integrals_1e<defs::ham_t, defs::isym_1e> ints1_t;
    typedef Integrals_2e<defs::ham_t, defs::isym_2e> ints2_t;
    ints1_t m_int_1;
    ints2_t m_int_2;

    ham_data::TermContribs m_contribs_1100;
    ham_data::TermContribs m_contribs_2200;
    ham_data::KramersAttributes m_kramers_attrs;

    AbinitioHamiltonian(size_t nelec, size_t nsite, bool spin_resolved, int ms2_restrict,
                        bool complex_valued, defs::inds site_irreps = {});

    AbinitioHamiltonian(const FcidumpHeader& header, bool spin_major, int ms2_restrict, int charge = 0);

    AbinitioHamiltonian(std::string fname, bool spin_major, int charge = 0):
            AbinitioHamiltonian(FcidumpHeader(fname), spin_major, charge){}

    explicit AbinitioHamiltonian(const fciqmc_config::FermionHamiltonian &opts) :
            AbinitioHamiltonian(opts.m_fcidump.m_path, opts.m_fcidump.m_spin_major, opts.m_charge) {}

    defs::ham_t get_coeff_1100(const size_t &i, const size_t &j) const override;

    defs::ham_t get_coeff_2200(const size_t &i, const size_t &j, const size_t &k, const size_t &l) const override;


    defs::ham_t get_element_0000(const field::FrmOnv &onv) const override;

    defs::ham_t get_element_1100(const field::FrmOnv &onv, const conn::FrmOnv &conn) const override {
        DEBUG_ASSERT_EQ(conn.exsig(), exsig_utils::ex_single, "expected 1100 (aka fermion single) exsig");
        const auto &ann = conn.m_ann[0];
        const auto &cre = conn.m_cre[0];

        defs::ham_t element = m_int_1(cre, ann);
        auto fn = [&](const size_t &ibit) {
            if (ibit != ann) element += m_int_2.phys_antisym_element(cre, ibit, ann, ibit);
        };
        onv.foreach(fn);
        return conn.phase(onv) ? -element : element;
    }

    defs::ham_t get_element_2200(const field::FrmOnv &onv, const conn::FrmOnv &conn) const override {
        DEBUG_ASSERT_EQ(conn.exsig(), exsig_utils::ex_double, "expected 2200 (aka fermion double) exsig");
        const auto element = m_int_2.phys_antisym_element(conn.m_cre[0], conn.m_cre[1], conn.m_ann[0], conn.m_ann[1]);
        return conn.phase(onv) ? -element : element;
    }

    buffered::FrmOnv guess_reference(const int &spin_level) const;

private:
    /**
     * if electrons are disabled, the m_nsite members of the integrals will be set to zero so as
     * not to allocate unnecessary memory
     * @return
     * number of distinctly specified (spin-)orbitals in the Hamiltonian defintion
     */
    size_t norb_distinct() const {
        return (1ul+m_int_1.m_spin_res)*m_nsite;
    }
};


#endif //M7_ABINITIOHAMILTONIAN_H
