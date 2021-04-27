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
 *
 *
 * Scalars associated with terms are called coefficients, and those associated with connections are called elements.
 *
 */
class FermionHamiltonian {
protected:
    const size_t m_nelec;
    const size_t m_nsite;
    const bool m_spin_conserving_1e, m_spin_conserving_2e;
    const bool m_complex_valued;
    const size_t m_int_2e_rank;
    defs::ham_t m_int_0;
    typedef Integrals_1e<defs::ham_t, defs::isym_1e> ints1_t;
    typedef Integrals_2e<defs::ham_t, defs::isym_2e> ints2_t;
    ints1_t m_int_1;
    ints2_t m_int_2;

public:
    FermionHamiltonian(const size_t &nelec, const size_t &nsite, bool spin_conserving_1e, bool spin_conserving_2e,
                       bool complex_valued, bool spin_resolved, size_t int_2e_rank);

    FermionHamiltonian(const FcidumpFileReader &file_reader);

    FermionHamiltonian(std::string fname, bool spin_major);

    FermionHamiltonian(const Options& opts): FermionHamiltonian(opts.fcidump_path, opts.fcidump_spin_major) {
        std::cout << opts.fcidump_path << std::endl;
    }

    defs::ham_comp_t get_energy(const fields::Onv<0> &fonv) const;

    defs::ham_t get_element_0(const defs::inds &occs, const size_t &nocc) const;

    defs::ham_t get_element_0(const OccupiedOrbitals &occs) const;

    defs::ham_t get_element_0(const fields::Onv<0> &fonv) const;

    defs::ham_t get_element_0(const conn::Antisym<0> &connection) const;

    defs::ham_t get_element_1(const conn::Antisym<0> &connection) const;

    defs::ham_t get_element_2(const size_t &i, const size_t &j, const size_t &k, const size_t &l) const;

    defs::ham_t get_element_2(const conn::Basic<0> &connection) const;

    defs::ham_t get_element_2(const conn::Antisym<0> &connection) const;

    defs::ham_t get_element(const conn::Antisym<0> &connection) const;

    defs::ham_t get_element(const fields::Onv<0> &bra, const fields::Onv<0> &ket) const;

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

    const size_t& int_2e_rank() const {
        return m_int_2e_rank;
    }

    buffered::FermionOnv guess_reference(const int &spin_level) const;
//
//    elements::FermionOnv refine_guess_reference(const views::FermionOnv &ref) const;
//
//    elements::FermionOnv choose_reference(const int &spin_level) const;
//
//    class DeterminantList : public MappedList<DeterminantElement> {
//    public:
//        DeterminantField determinant;
//
//        DeterminantList(std::string name, size_t nsite, size_t nbucket) :
//                MappedList(name, determinant, nbucket),
//                determinant(this, 1, nsite){}
//    };
//
//    void generate_ci_space(WalkerTable* list, RankAllocator<DeterminantElement>& ra, const int &spin_level) const;

};

#endif //M7_FERMIONHAMILTONIAN_H
