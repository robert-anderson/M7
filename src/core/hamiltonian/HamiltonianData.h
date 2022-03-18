//
// Created by rja on 21/08/2021.
//

#ifndef M7_HAMILTONIANDATA_H
#define M7_HAMILTONIANDATA_H


#include <util/utils.h>
#include <parallel/MPIAssert.h>


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

        void set_nonzero(size_t exsig);

        bool is_nonzero(size_t exsig) const;

    };

    struct FrmModelAttributes {

        bool m_nn_only_singles = true;
        /**
         * same as above, but with (1D) periodic boundary conditions
         */
        bool m_nnp_only_singles = true;

        /**
         * are (2, 2)-rank contribs due to (0, 0)-excitations nonzero only when all operator indices refer to the same site?
         * assume this is the case unless counterexample found in loop over nonzero elements
         */
        bool m_on_site_only_doubles = true;

        bool is_hubbard_1d() const;

        bool is_hubbard_1d_pbc() const;
        /**
         * @param iorb
         *  orbital index (site if integrals spin resolved, else spin oribtal index)
         * @return
         *  site index
         */
        static size_t iorb_to_isite(const size_t &iorb, const size_t &nsite);

        /**
         * @param iorb
         * @param jorb
         * @param periodic
         *  if true, lowest and highest site indices count as neighbors.
         * @return
         * true if orbitals are nearest neighbors
         */
        static bool nearest_neighbors(const size_t &nsite, const size_t &iorb, const size_t &jorb, bool periodic);

        bool on_site(const size_t &nsite, const size_t &iorb, const size_t &jorb, const size_t &korb, const size_t &lorb) const;

        void nonzero(const size_t & nsite, const size_t & i, const size_t & j);

        /**
         * non-zero integral reached with chemical notation (ij|kl)
         */
        void nonzero(const size_t & nsite, const size_t & i, const size_t & j, const size_t & k, const size_t & l);
    };

    struct KramersAttributes {
        bool m_conserving_singles = true;
        bool m_conserving_double = true;

        bool conserving() const;
    };
}

#endif //M7_HAMILTONIANDATA_H
