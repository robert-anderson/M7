//
// Created by rja on 21/08/2021.
//

#ifndef M7_HAMILTONIANDATA_H
#define M7_HAMILTONIANDATA_H


#include <M7_lib/util/utils.h>
#include <M7_lib/parallel/MPIAssert.h>


namespace ham_data {
    using namespace exsig_utils;

    class TermContribs {

        const size_t m_ranksig;
        const size_t m_basesig;
        const size_t m_nexsig_contrib_frm, m_nexsig_contrib_bos;

        std::vector<bool> m_exsig_nonzero;

        size_t ind(size_t exsig) const;


    public:
        TermContribs(size_t ranksig);

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

        void set_nonzero(size_t exsig);

        bool is_nonzero(size_t exsig) const;

    };

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
}

#endif //M7_HAMILTONIANDATA_H
