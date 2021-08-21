//
// Created by rja on 27/02/2020.
//

#ifndef M7_FERMIONHAMILTONIAN_H
#define M7_FERMIONHAMILTONIAN_H

#include <cstddef>
#include <src/core/basis/DecodedDeterminant.h>
#include <src/core/connection/Connections.h>
#include <src/core/io/Options.h>
#include <src/core/config/FciqmcConfig.h>
#include "src/core/integrals/Integrals_1e.h"
#include "src/core/integrals/Integrals_2e.h"
#include "src/core/table/BufferedFields.h"
#include "HamiltonianData.h"

/**
 * All interactions between the fermionic parts of ONVs are described in this class.
 */
struct FermionHamiltonian {

    const size_t m_nelec;
    const size_t m_nsite;
    const bool m_spin_conserving_1e, m_spin_conserving_2e;
    const bool m_complex_valued;
    const size_t m_int_2e_rank;

public:
    defs::ham_t m_int_0;
    typedef Integrals_1e<defs::ham_t, defs::isym_1e> ints1_t;
    typedef Integrals_2e<defs::ham_t, defs::isym_2e> ints2_t;
    ints1_t m_int_1;
    ints2_t m_int_2;

    /**
     * stores whether connected elements are nonzero based on contrib case
     */
    std::array<bool, ham_data::contrib_count> m_nonzero_contribs;
    /**
     * are (1, 1)-rank contribs due to (1, 1)-excitations nonzero only for (1D) nearest neighbors?
     * assume this is the case unless counterexample found in loop over nonzero elements
     */
    bool m_nn_only_1111 = true;
    /**
     * same as above, but with (1D) periodic boundary conditions
     */
    bool m_nnp_only_1111 = true;
    /**
     * are (2, 2)-rank contribs due to (0, 0)-excitations nonzero only when all operator indices refer to the same site?
     * assume this is the case unless counterexample found in loop over nonzero elements
     */
    bool m_on_site_only_0022 = true;

public:

    /**
     * @param iorb
     *  orbital index (site if integrals spin resolved, else spin oribtal index)
     * @return
     *  site index
     */
    size_t iorb_to_isite(const size_t &iorb) const {
        return iorb < m_nsite ? iorb : iorb + m_nsite;
    }

    /**
     * @param iorb
     * @param jorb
     * @param periodic
     *  if true, lowest and highest site indices count as neighbors.
     * @return
     * true if orbitals are nearest neighbors
     */
    bool nearest_neighbors(const size_t &iorb, const size_t &jorb, bool periodic) const {
        auto isite = iorb_to_isite(iorb);
        auto jsite = iorb_to_isite(jorb);
        if (isite + 1 == jsite || jsite + 1 == isite) return true;
        if (periodic) {
            return (isite == 0 && jsite == m_nsite - 1) || (isite == m_nsite - 1 && jsite == 0);
        }
        return false;
    }

    bool on_site(const size_t &iorb, const size_t &jorb, const size_t &korb, const size_t &lorb) const {
        auto isite = iorb_to_isite(iorb);
        auto jsite = iorb_to_isite(jorb);
        auto ksite = iorb_to_isite(korb);
        auto lsite = iorb_to_isite(lorb);
        return isite == jsite && jsite == ksite && ksite == lsite;
    }

    FermionHamiltonian(const size_t &nelec, const size_t &nsite, bool spin_conserving_1e, bool spin_conserving_2e,
                       bool complex_valued, bool spin_resolved, size_t int_2e_rank);

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
        DEBUG_ASSERT_EQ(conn.exsig(), conn_utils::encode_exsig(1,1,0,0), "expected 1100 exsig");
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
        DEBUG_ASSERT_EQ(conn.exsig(), conn_utils::encode_exsig(2,2,0,0), "expected 2200 exsig");
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

public:
    bool is_hubbard_1d() const {
        return m_on_site_only_0022 && m_nn_only_1111;
    }

    bool is_hubbard_1d_pbc() const {
        return m_on_site_only_0022 && m_nnp_only_1111;
    }

    bool spin_conserving() const {
        return m_spin_conserving_1e && m_spin_conserving_2e;
    }

    buffered::FrmOnv guess_reference(const int &spin_level) const;

    /**
     * set the referenced ONV object to the assumed Hartree--Fock determinant within the given spin sector
     * @param onv
     *  target onv object
     * @param spin
     *  spin (MS) number
     */
    void set_hf_mbf(field::FrmOnv &onv, int spin) const {
        auto nalpha = ci_utils::nalpha(m_nelec, spin);
        auto nbeta = ci_utils::nbeta(m_nelec, spin);
        DEBUG_ASSERT_EQ(nalpha + nbeta, m_nelec, "inconsistent na, nb, nelec");
        onv.zero();
        for (size_t i = 0ul; i < nalpha; ++i) onv.set({0, i});
        for (size_t i = 0ul; i < nbeta; ++i) onv.set({1, i});
    }

};

#endif //M7_FERMIONHAMILTONIAN_H
