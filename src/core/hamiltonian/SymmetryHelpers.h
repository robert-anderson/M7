//
// Created by rja on 24/05/2021.
//

#ifndef M7_SYMMETRYHELPERS_H
#define M7_SYMMETRYHELPERS_H

#include <src/core/basis/Connections.h>
#include <src/core/basis/DecodedDeterminant.h>
#include <src/core/table/BufferedFields.h>

/**
 * symmetry helpers use symmetries of the hamiltonian to reduce redundant looping in:
 *  1. the computation of hamiltonian matrix elements
 *  2. deterministic enumeration of connections to a many-body basis function
 *
 * 1. is handled by the classes in the hams namespace by default, but for some lattice models it may be more efficient
 *    to compute these in such a way that the runtime-detected symmetry of the Hamiltonian is respected.
 * 2. is only handled here. it is useful for exact propagation, test code and sign problem free energy estimators, and
 *    will likely have other applications in future
 */


struct FermionHamiltonian;
struct FermiBosHamiltonian;

namespace ham_sym_helpers {

    /**
     * All interactions between the fermionic parts of ONVs are described in this namespace.
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
     */
    // no-sym base class
    struct Fermion {

        typedef std::function<void(const conn::Antisym<0> &, const fields::Onv<0> &, const defs::ham_t &)> body_fn_t;

        const FermionHamiltonian &m_ham;
        const size_t m_nsite;

        mutable defs::ham_t m_helement_work;
        mutable conn::Antisym<0> m_conn_work;
        mutable OccupiedOrbitals m_occ_work;
        mutable VacantOrbitals m_vac_work;
        mutable buffered::Onv<0> m_onv_work;

        Fermion(const FermionHamiltonian &ham);

        virtual ~Fermion(){};

        virtual defs::ham_t get_element(const conn::Antisym<0> &conn) const;

        defs::ham_comp_t get_energy(const fields::Onv<0> &onv) const;

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

        void perform_diagonal(const fields::Onv<0> &src_onv,
                              const body_fn_t &body, bool get_h, bool h_nonzero_only) const;

        void perform_single(const fields::Onv<0> &src_onv, const size_t& occ, const size_t& vac,
                            const body_fn_t &body, bool get_h, bool h_nonzero_only) const;

        void perform_double(const fields::Onv<0> &src_onv,
                            const size_t& occ1, const size_t& occ2,
                            const size_t& vac1, const size_t& vac2,
                            const body_fn_t &body, bool get_h, bool h_nonzero_only) const;

        void foreach_connection_singles(const fields::Onv<0> &src_onv,
                                        const defs::inds &occs, const defs::inds &vacs,
                                        const body_fn_t &body, bool get_h, bool h_nonzero_only) const;

        void foreach_connection_subset_doubles(const fields::Onv<0> &src_onv,
                                               const defs::inds &occs1, const defs::inds &occs2,
                                               const defs::inds &vacs1, const defs::inds &vacs2,
                                               const body_fn_t &body, bool get_h, bool h_nonzero_only) const;

        virtual void foreach_connection_subset_doubles(const fields::Onv<0> &src_onv,
                                                       const defs::inds &occs, const defs::inds &vacs,
                                                       const body_fn_t &body, bool get_h, bool h_nonzero_only) const;

        void foreach_connection_subset(const fields::Onv<0> &src_onv,
                                       const defs::inds &occs1, const defs::inds &occs2,
                                       const defs::inds &vacs1, const defs::inds &vacs2,
                                       const body_fn_t &body, bool get_h, bool h_nonzero_only) const;

        virtual void foreach_connection_subset(const fields::Onv<0> &src_onv,
                                               const defs::inds &occs, const defs::inds &vacs,
                                               const body_fn_t &body, bool get_h, bool h_nonzero_only) const;

    public:
        virtual void foreach_connection(const fields::Onv<0> &src_onv, const body_fn_t &body,
                                        bool get_h, bool h_nonzero_only, bool include_diagonal) const;
    };


    // no-sym base class
    struct FermiBos {

        typedef std::function<void(const conn::Antisym<1> &, const fields::Onv<1> &, const defs::ham_t &)> body_fn_t;

        const FermiBosHamiltonian &m_ham;
        const size_t m_nsite;

        mutable defs::ham_t m_helement_work;
        mutable conn::Antisym<1> m_conn_work;
        mutable OccupiedOrbitals m_occ_work;
        mutable VacantOrbitals m_vac_work;
        mutable buffered::Onv<1> m_onv_work;

        FermiBos(const FermiBosHamiltonian &ham);

        virtual ~FermiBos(){};

    private:
        bool update_helement(bool get_h, bool h_nonzero_only) const {
            if (get_h || h_nonzero_only) m_helement_work = get_element(m_conn_work);
            else m_helement_work = 0.0;
            if (h_nonzero_only) return !consts::float_is_zero(m_helement_work);
            return true;
        }

    public:

        virtual defs::ham_t get_element(const conn::Antisym<1> &conn) const;

        defs::ham_comp_t get_energy(const fields::Onv<1> &onv) const;

        virtual void foreach_connection(const fields::Onv<1> &src_onv, const body_fn_t &body,
                                        bool get_h, bool h_nonzero_only, bool include_diagonal) const;

    };
}


#endif //M7_SYMMETRYHELPERS_H
