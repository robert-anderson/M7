//
// Created by rja on 27/02/2020.
//

#ifndef M7_FERMIONHAMILTONIAN_H
#define M7_FERMIONHAMILTONIAN_H

#include <cstddef>
#include <src/core/basis/FermionOnvConnection.h>
#include <src/core/basis/DecodedDeterminant.h>
#include <src/core/basis/Connections.h>
#include <src/core/io/Options.h>
#include "src/core/integrals/Integrals_1e.h"
#include "src/core/integrals/Integrals_2e.h"
#include "src/core/table/BufferedFields.h"
#include "HamiltonianData.h"
#include "SymmetryHelpers.h"


/**
 * All interactions between the fermionic parts of ONVs are described in this class.
 *
 * The rank of a normal-ordered fermion operator product is denoted by two integers as (nc, na) where nc is the number
 * of creation operators in the product (which must be distinct by PEP) and na is the number of annihilation operators
 * in the product with the same restriction.
 *
 * The standard quantum-chemical Hamiltonian has terms of three (fermion number conserving) ranks
 * (0, 0), (1, 1), (2, 2)
 *
 * However, higher values of nc and na are conceivably of interest in model studies, and although they are not included
 * in the current implementation, they can be easily incorporated with the definition of new methods following the
 * naming convention set out here.
 *
 * On the other hand, fermion number non-conserving interactions are implemented, which have ranks denoted by
 * (0, 1), (1, 2), (0, 2),  (1, 0), (2, 1), (2, 0)
 * effecting particle number sector changes (from ket to bra) +1, +1, +2, -1, -1, -2 respectively.
 *
 * Terms of a given rank are subdivided into multiple excitation levels, and each ONV connection is identified with one
 * of these. The notation here is similar to that of the H terms:
 *
 * Excitation level of connection  |  May require H element contributions from terms of rank:
 * -----------------------------------------------------------------------------
 *      (0, 0)                     |  (0, 0), (1, 1), (2, 2)
 *      (1, 1)                     |  (1, 1), (2, 2)
 *      (2, 2)                     |  (2, 2)
 *      (0, 1)                     |  (0, 1), (1, 2)
 *      (1, 2)                     |  (1, 2)
 *      (0, 2)                     |  (0, 2)
 *      (1, 0)                     |  (1, 0), (2, 1)
 *      (2, 1)                     |  (2, 1)
 *      (2, 0)                     |  (2, 0)
 * ------------------------------------------------------------------------------
 *
 * Clearly, only four of the presently handled excitation levels require contributions from higher-rank terms in the
 * Hamiltonian. Some models may possess symmetry properties which eliminate any or all of the contributions from an
 * excitation level
 * e.g. N-dimensional Hubbard model has
 *  * no rank (1, 1) contribution to exlvl (0, 0) - (no "hopping" to same spin orbital)
 *  * no rank (2, 2) contribution to exlvl (1, 1) - (two-body interaction is diagonal)
 *  * no exlvl (2, 2) contributions at all - (two-body interaction is diagonal)
 *
 * A fermion H element contribution is given a 4-digit identifier wxyz,
 * where (w, x) is the excitation level, and (y, z) is the term rank.
 *
 *
 * Scalars associated with terms are called coefficients, and those associated with connections are called elements.
 *
 * methods to return term coefficients for every implemented rank (get_coeff_)
 * methods to compute matrix elements (get_element_)
 *
 *
 */
struct FermionHamiltonian {

protected:
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
    const std::unique_ptr<ham_sym_helpers::Fermion> m_sym_helper;

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
            return (isite == 0 && jsite == m_nsite-1) || (isite == m_nsite-1 && jsite == 0);
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

    FermionHamiltonian(const Options &opts) : FermionHamiltonian(opts.fcidump_path, opts.fcidump_spin_major) {}

    defs::ham_t get_element_0(const defs::inds &occs, const size_t &nocc) const;

    defs::ham_t get_element_0(const OccupiedOrbitals &occs) const;

    defs::ham_t get_element_0(const fields::Onv<0> &fonv) const;

    defs::ham_t get_element_0(const conn::Antisym<0> &connection) const;

    defs::ham_t get_element_1(const conn::Antisym<0> &connection) const;

    defs::ham_t get_element_2(const size_t &i, const size_t &j, const size_t &k, const size_t &l) const;

    defs::ham_t get_element_2(const conn::Basic<0> &connection) const;

    defs::ham_t get_element_2(const conn::Antisym<0> &connection) const;

    defs::ham_t get_element(const fields::Onv<0> &bra, const fields::Onv<0> &ket) const;

private:
    defs::ham_t get_element_tag(const conn::Antisym<0> &conn, dispatch_utils::BoolTag<false> sym_opts) const {
        switch (conn.nexcit()) {
            case 0:
                return get_element_0(conn);
            case 1: ASSERT(conn.ncom() + conn.nexcit() == nelec());
                return get_element_1(conn);
            case 2:
                return get_element_2(conn);
            default:
                return 0;
        }
    }
    defs::ham_t get_element_tag(const conn::Antisym<0> &conn, dispatch_utils::BoolTag<true> sym_opts) const {
        return m_sym_helper->get_element(conn);
    }

    defs::ham_comp_t get_energy_tag(const fields::Onv<0> &onv, dispatch_utils::BoolTag<false> sym_opts) const {
        return consts::real(get_element_0(onv));
    }
    defs::ham_comp_t get_energy_tag(const fields::Onv<0> &onv, dispatch_utils::BoolTag<true> sym_opts) const {
        return m_sym_helper->get_energy(onv);
    }

public:
    defs::ham_t get_element(const conn::Antisym<0> &conn) const {
        return get_element_tag(conn, dispatch_utils::BoolTag<defs::enable_optim_for_lattice_ham>());
    }

    defs::ham_comp_t get_energy(const fields::Onv<0> &onv) const {
        return get_energy_tag(onv, dispatch_utils::BoolTag<defs::enable_optim_for_lattice_ham>());
    }

    bool is_hubbard() const {
        return m_on_site_only_0022 && m_nn_only_1111;
    }

    bool is_hubbard_pbc() const {
        return m_on_site_only_0022 && m_nnp_only_1111;
    }

    size_t nci() const {
        return ci_utils::fermion_dim(nsite(), nelec());
    }

    const size_t &nsite() const {
        return m_nsite;
    }

    bool spin_conserving_1e() const {
        return m_spin_conserving_1e;
    }

    bool spin_conserving_2e() const {
        return m_spin_conserving_2e;
    }

    bool spin_conserving() const {
        return m_spin_conserving_1e && m_spin_conserving_2e;
    }

    const size_t &nelec() const {
        return m_nelec;
    }

    const bool &complex_valued() const {
        return m_complex_valued;
    }

    const size_t &int_2e_rank() const {
        return m_int_2e_rank;
    }

    buffered::FermionOnv guess_reference(const int &spin_level) const;


    void foreach_connection(const fields::Onv<0> &src_onv, const ham_sym_helpers::Fermion::body_fn_t &body,
                            bool get_h, bool h_nonzero_only, bool include_diagonal) const {
        m_sym_helper->foreach_connection(src_onv, body, get_h, h_nonzero_only, include_diagonal);
    }

    /**
     * set the referenced ONV object to the assumed Hartree--Fock determinant within the given spin sector
     * @param onv
     *  target onv object
     * @param spin
     *  spin (MS) number
     */
    void set_hf_onv(fields::FermionOnv& onv, int spin) const {
        auto nalpha = ci_utils::nalpha(nelec(), spin);
        auto nbeta = ci_utils::nbeta(nelec(), spin);
        MPI_ASSERT(nalpha+nbeta==nelec(), "inconsistent na, nb, nelec");
        onv.zero();
        for (size_t i=0ul; i<nalpha; ++i) onv.set({0, i});
        for (size_t i=0ul; i<nbeta; ++i) onv.set({1, i});
    }

};

#endif //M7_FERMIONHAMILTONIAN_H
