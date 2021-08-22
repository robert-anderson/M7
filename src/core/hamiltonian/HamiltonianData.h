//
// Created by rja on 21/08/2021.
//

#ifndef M7_HAMILTONIANDATA_H
#define M7_HAMILTONIANDATA_H


#include <src/core/util/utils.h>
#include <src/core/parallel/MPIAssert.h>


namespace ham_data {
    using namespace conn_utils;

    class TermContribs {

        const size_t m_ranksig;
        const size_t m_basesig;
        const size_t m_nexsig_contrib_frm, m_nexsig_contrib_bos;

        std::vector<bool> m_exsig_nonzero;

        size_t ind(size_t exsig) const {
            if (!contribs_to(exsig, m_ranksig)) return ~0ul;
            size_t ifrm = decode_nfrm_cre(exsig);
            DEBUG_ASSERT_LT(ifrm, m_nexsig_contrib_frm, "invalid number of like-indexed fermion operators");
            size_t ibos = decode_nbos_cre(exsig);
            DEBUG_ASSERT_LT(ibos, m_nexsig_contrib_bos, "invalid number of like-indexed boson operators");
            return ifrm*m_nexsig_contrib_bos+ibos;
        }


    public:
        TermContribs(size_t ranksig): m_ranksig(ranksig), m_basesig(base_exsig(ranksig)),
                                      m_nexsig_contrib_frm(ncontrib_frm(ranksig)), m_nexsig_contrib_bos(ncontrib_bos(ranksig)),
                                      m_exsig_nonzero(m_nexsig_contrib_frm*m_nexsig_contrib_bos, false){}

        void set_nonzero(size_t exsig) {
            m_exsig_nonzero[ind(exsig)] = true;
        }

        bool is_nonzero(size_t exsig) const {
            return m_exsig_nonzero[ind(exsig)];
        }

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

        bool is_hubbard_1d() const {
            return m_on_site_only_doubles && m_nn_only_singles;
        }

        bool is_hubbard_1d_pbc() const {
            return m_on_site_only_doubles && m_nnp_only_singles;
        }
        /**
         * @param iorb
         *  orbital index (site if integrals spin resolved, else spin oribtal index)
         * @return
         *  site index
         */
        static size_t iorb_to_isite(const size_t &iorb, const size_t &nsite) {
            return iorb < nsite ? iorb : iorb + nsite;
        }

        /**
         * @param iorb
         * @param jorb
         * @param periodic
         *  if true, lowest and highest site indices count as neighbors.
         * @return
         * true if orbitals are nearest neighbors
         */
        static bool nearest_neighbors(const size_t &nsite, const size_t &iorb, const size_t &jorb, bool periodic) {
            auto isite = iorb_to_isite(iorb, nsite);
            auto jsite = iorb_to_isite(jorb, nsite);
            if (isite + 1 == jsite || jsite + 1 == isite) return true;
            if (periodic) {
                return (isite == 0 && jsite == nsite - 1) || (isite == nsite - 1 && jsite == 0);
            }
            return false;
        }

        bool on_site(const size_t &nsite, const size_t &iorb, const size_t &jorb, const size_t &korb, const size_t &lorb) const {
            auto isite = iorb_to_isite(iorb, nsite);
            auto jsite = iorb_to_isite(jorb, nsite);
            auto ksite = iorb_to_isite(korb, nsite);
            auto lsite = iorb_to_isite(lorb, nsite);
            return isite == jsite && jsite == ksite && ksite == lsite;
        }

        void nonzero(const size_t & nsite, const size_t & i, const size_t & j){
            if (i==j){
                // one-electron term is not purely off-diagonal
                m_nn_only_singles = false;
                m_nnp_only_singles = false;
            }
            else {
                if (!nearest_neighbors(nsite, i, j, false)) m_nn_only_singles = false;
                if (!nearest_neighbors(nsite, i, j, true)) m_nnp_only_singles = false;
            }
        }

        /**
         * non-zero integral reached with chemical notation (ij|kl)
         */
        void nonzero(const size_t & nsite, const size_t & i, const size_t & j, const size_t & k, const size_t & l){
            if (!(i==j && j==k && k==l)) m_on_site_only_doubles = false;
        }
    };

    struct KramersAttributes {
        bool m_conserving_singles = true;
        bool m_conserving_double = true;

        bool conserving() const {
            return m_conserving_singles && m_conserving_double;
        }
    };
}

#endif //M7_HAMILTONIANDATA_H
