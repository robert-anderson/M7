//
// Created by rja on 27/02/2020.
//

#ifndef M7_FERMIONHAMILTONIAN_H
#define M7_FERMIONHAMILTONIAN_H

#include <cstddef>
#include <src/core/basis/DecodedDeterminants.h>
#include <src/core/connection/Connections.h>
#include <src/core/io/Options.h>
#include <src/core/config/FciqmcConfig.h>
#include "src/core/integrals/Integrals_1e.h"
#include "src/core/integrals/Integrals_2e.h"
#include "src/core/table/BufferedFields.h"
#include "HamiltonianData.h"


/**
 * All interactions between the fermionic parts of MBFs are described in this class.
 */
struct FermionHamiltonian {

    const size_t m_nelec;
    const size_t m_nsite;
    const bool m_complex_valued;
    AbelianGroupMap m_point_group_map;

    defs::ham_t m_int_0 = 0.0;
    typedef Integrals_1e<defs::ham_t, defs::isym_1e> ints1_t;
    typedef Integrals_2e<defs::ham_t, defs::isym_2e> ints2_t;
    ints1_t m_int_1;
    ints2_t m_int_2;

    ham_data::TermContribs m_contribs_1100;
    ham_data::TermContribs m_contribs_2200;
    ham_data::FrmModelAttributes m_model_attrs;
    ham_data::KramersAttributes m_kramers_attrs;

    FermionHamiltonian(size_t nelec, size_t nsite, bool complex_valued,
                       bool spin_resolved, defs::inds site_irreps={});

    FermionHamiltonian(const FcidumpFileReader &file_reader);

    FermionHamiltonian(std::string fname, bool spin_major);

    FermionHamiltonian(const fciqmc_config::Hamiltonian &opts) :
            FermionHamiltonian(opts.m_fcidump.m_path, opts.m_fcidump.m_spin_major) {}

    defs::ham_t get_element_0000(const field::FrmOnv &onv) const;

    defs::ham_t get_element(const field::FrmOnv &onv) const {
        return get_element_0000(onv);
    }

    defs::ham_comp_t get_energy(const field::FrmOnv &onv) const {
        return consts::real(get_element_0000(onv));
    }

    defs::ham_t get_element_1100(const field::FrmOnv &onv, const conn::FrmOnv &conn) const {
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

    defs::ham_t get_element_2200(const size_t &i, const size_t &j, const size_t &k, const size_t &l) const;

    defs::ham_t get_element_2200(const field::FrmOnv &onv, const conn::FrmOnv &conn) const {
        DEBUG_ASSERT_EQ(conn.exsig(), exsig_utils::ex_double, "expected 2200 (aka fermion double) exsig");
        const auto element = get_element_2200(conn.m_cre[0], conn.m_cre[1], conn.m_ann[0], conn.m_ann[1]);
        return conn.phase(onv) ? -element : element;
    }

    defs::ham_t get_element(const field::FrmOnv &onv, const conn::FrmOnv &conn) const {
        switch (conn.size()) {
            case 0: return get_element_0000(onv);
            case 2: return get_element_1100(onv, conn);
            case 4: return get_element_2200(onv, conn);
            default: return 0.0;
        }
    }

    size_t nci() const {
        return ci_utils::fermion_dim(m_nsite, m_nelec);
    }

    buffered::FrmOnv guess_reference(const int &spin_level) const;

    /**
     * set the referenced ONV object to the assumed Hartree--Fock determinant within the given spin sector
     * @param onv
     *  target onv object
     * @param spin
     *  spin (MS) number
     */
    void set_hf_mbf(field::FrmOnv &onv, int spin) const;

    /**
     * output some useful logs identifying the kind of H detected
     */
    void log_data() const;

    bool is_hubbard_1d() const {
        return m_model_attrs.is_hubbard_1d();
    }
};

#endif //M7_FERMIONHAMILTONIAN_H
