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

    /**
     * Base class virtual methods take no symmetry into account and is fermion number conserving
     */
    struct Terms {
        typedef std::function<void(const conn::Antisym<0> &, const fields::FermionOnv &,
                                   const defs::ham_t &)> body_fn_t;
        const FermionHamiltonian &m_ham;

        mutable defs::ham_t m_helement_work;
        mutable conn::Antisym<0> m_conn_work;
        mutable OccupiedOrbitals m_occ_work;
        mutable VacantOrbitals m_vac_work;
        mutable buffered::FermionOnv m_onv_work;

        Terms(const FermionHamiltonian &ham) : m_ham(ham), m_conn_work(ham.nsite()),
                                               m_occ_work(ham.nsite()), m_vac_work(ham.nsite()),
                                               m_onv_work(ham.nsite()) {}

        virtual defs::ham_t get_element_00(const defs::inds &occs, const size_t &nocc) const;

        virtual defs::ham_t get_element_11(const conn::Antisym<0> &conn) const;

        virtual defs::ham_t get_element_22(const conn::Antisym<0> &conn) const;

        virtual defs::ham_t get_element_01(const conn::Antisym<0> &conn) const;

        virtual defs::ham_t get_element_12(const conn::Antisym<0> &conn) const;

        virtual defs::ham_t get_element_02(const conn::Antisym<0> &conn) const;

        virtual defs::ham_t get_element_10(const conn::Antisym<0> &conn) const;

        virtual defs::ham_t get_element_21(const conn::Antisym<0> &conn) const;

        virtual defs::ham_t get_element_20(const conn::Antisym<0> &conn) const;

        defs::ham_t get_element(const conn::Antisym<0> &conn) const;

    private:
        /**
         * decide whether we need to compute full H element, and if so, whether the should the loop body be called
         * @param get_h
         *  if true, the calling scope wants access to the H matrix element
         * @param h_nonzero_only
         *  if true, the calling scope only wants the loop body to be called when the matrix element is nonzero
         * @return
         *  true if the loop body should be called
         */
        bool update_helement(bool get_h, bool h_nonzero_only) const;

    protected:

        void perform_single(const fields::FermionOnv &src_onv, const size_t& occ, const size_t& vac,
                    const body_fn_t &body, bool get_h, bool h_nonzero_only) const;

        void perform_double(const fields::FermionOnv &src_onv,
                    const size_t& occ1, const size_t& occ2,
                    const size_t& vac1, const size_t& vac2,
                    const body_fn_t &body, bool get_h, bool h_nonzero_only) const;

        void foreach_connection_singles(const fields::FermionOnv &src_onv,
                                       const defs::inds &occs, const defs::inds &vacs,
                                       const body_fn_t &body, bool get_h, bool h_nonzero_only) const;

        void foreach_connection_subset_doubles(const fields::FermionOnv &src_onv,
                                       const defs::inds &occs1, const defs::inds &occs2,
                                       const defs::inds &vacs1, const defs::inds &vacs2,
                                       const body_fn_t &body, bool get_h, bool h_nonzero_only) const;

        virtual void foreach_connection_subset_doubles(const fields::FermionOnv &src_onv,
                                               const defs::inds &occs, const defs::inds &vacs,
                                               const body_fn_t &body, bool get_h, bool h_nonzero_only) const;

        void foreach_connection_subset(const fields::FermionOnv &src_onv,
                                               const defs::inds &occs1, const defs::inds &occs2,
                                               const defs::inds &vacs1, const defs::inds &vacs2,
                                               const body_fn_t &body, bool get_h, bool h_nonzero_only) const;

        virtual void foreach_connection_subset(const fields::FermionOnv &src_onv,
                                               const defs::inds &occs, const defs::inds &vacs,
                                               const body_fn_t &body, bool get_h, bool h_nonzero_only) const;

    public:
        /**
         * A loop over ONVs connected to src_onv
         * @param src_onv
         * @param body
         * @param get_h
         * @param h_nonzero_only
         */
        virtual void foreach_connection(const fields::FermionOnv &src_onv, const body_fn_t &body,
                                        bool get_h, bool h_nonzero_only) const;
    };

    struct SpinTerms : Terms {
        mutable SpinOccupiedOrbitals m_spin_occ_work;
        mutable SpinVacantOrbitals m_spin_vac_work;

        SpinTerms(const FermionHamiltonian &ham);

        virtual void foreach_connection(const fields::FermionOnv &src_onv, const body_fn_t &body,
                                        bool get_h, bool h_nonzero_only) const;
    };


    struct Hubbard1DTerms : SpinTerms {

        Hubbard1DTerms(const FermionHamiltonian &ham);

        virtual void foreach_connection(const fields::FermionOnv &src_onv, const body_fn_t &body,
                                        bool get_h, bool h_nonzero_only) const;
    };

    struct Hubbard1DPbcTerms : SpinTerms {

        virtual void foreach_connection(const fields::FermionOnv &src_onv, const body_fn_t &body,
                                        bool get_h, bool h_nonzero_only) const;
    };

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
    const std::unique_ptr<Terms> m_terms;

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
            return (isite == 0 && jsite == m_nsite) || (isite == m_nsite && jsite == 0);
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

//    FermionHamiltonianTerms decide_terms() const {
//        return FermionHamiltonianTerms(*this);
//    }

    FermionHamiltonian(const size_t &nelec, const size_t &nsite, bool spin_conserving_1e, bool spin_conserving_2e,
                       bool complex_valued, bool spin_resolved, size_t int_2e_rank);

    FermionHamiltonian(const FcidumpFileReader &file_reader);

    FermionHamiltonian(std::string
                       fname,
                       bool spin_major
    );

    FermionHamiltonian(const Options &opts) : FermionHamiltonian(opts.fcidump_path, opts.fcidump_spin_major) {
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

    const size_t &int_2e_rank() const {
        return m_int_2e_rank;
    }

    buffered::FermionOnv guess_reference(const int &spin_level) const;


    void foreach_connection(const fields::FermionOnv &src_onv, const Terms::body_fn_t &body,
                            bool get_h, bool h_nonzero_only) const {
        m_terms->foreach_connection(src_onv, body, get_h, h_nonzero_only);
    }

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
