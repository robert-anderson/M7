//
// Created by Robert J. Anderson on 21/08/2021.
//

#ifndef M7_HAMILTONIANDATA_H
#define M7_HAMILTONIANDATA_H


#include "M7_lib/util/Exsig.h"
#include <M7_lib/parallel/MPIAssert.h>


namespace ham {
    using namespace exsig;

    class TermContribs {

        const uint_t m_ranksig;
        const uint_t m_basesig;
        const uint_t m_nexsig_contrib_frm, m_nexsig_contrib_bos;

        std::vector<bool> m_exsig_nonzero;

        uint_t ind(uint_t exsig) const;


    public:
        TermContribs(uint_t ranksig);

        TermContribs(const TermContribs& other);

        TermContribs& operator=(const TermContribs& other);

        /**
         * ctor to combine term contribs from summed hamiltonians
         * @param contribs_1
         *  zero/non-zero status of contributions from one hamiltonian
         * @param contribs_2
         *  zero/non-zero status of contributions from another hamiltonian
         */
        TermContribs(const TermContribs& contribs_1, const TermContribs& contribs_2);

        void set_nonzero(uint_t exsig);

        bool is_nonzero(uint_t exsig) const;

        bool any_nonzero() const;

    };

    /**
     * Scalar wavefunction theories are made compatible with the double degeneracy of electronic wavefunctions by
     * the ad hoc separation into two components with alpha and beta z-axis projection of spin, making the total Sz
     * trivially conserved. In relativistic theories however, terms expressed in terms of spin and spatial degrees of
     * freedom are not multiplicatively separable, and so Sz is not conserved. A notion of double degeneracy is
     * recovered by the Kramers time reversal symmetry operation, and in a Kramers-restricted basis (subject to certain
     * double point group symmetry considerations) the 4-component analogue of many body 2*Sz i.e. (Kramers+ - Kramers-)
     * can in some cases be conserved, but this is on a per-excitation signature basis.
     * see https://aip.scitation.org/doi/pdf/10.1063/5.0029863, Table 1 for a breakdown of these cases
     */
    struct KramersAttributes {
        bool m_conserving_singles = true;
        bool m_conserving_doubles = true;

        KramersAttributes(){}

        /**
         * ctor to combine kramers conservation attributes from summed hamiltonians
         * @param attrs_1
         *  kramers conservation/non-conservation status of one hamiltonian
         * @param attrs_2
         *  kramers conservation/non-conservation status of another hamiltonian
         */
        KramersAttributes(const KramersAttributes& attrs_1, const KramersAttributes& attrs_2);

        bool conserving() const;
    };

    static constexpr double element_tol() {return 1e-12;}
    /**
     * @param elem
     *  hamiltonian matrix element
     * @return
     *  true if the magnitude of elem exceeds hard-coded minimum
     */
    static bool is_significant(ham_t elem) {
        return fptol::nearly_zero(elem, element_tol());
    }
}

#endif //M7_HAMILTONIANDATA_H
