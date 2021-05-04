//
// Created by rja on 28/04/2021.
//

#ifndef M7_HAMILTONIANDATA_H
#define M7_HAMILTONIANDATA_H

#include "src/defs.h"

namespace ham_data {

    /**
     * keeping this number low allows the dispatch of rank-specific methods by rank-generic methods to be done via
     * switch statements and
     */
    constexpr size_t c_ndigit_rank_max = 2;
    constexpr size_t c_ndigit_rank_mask = (1ul<<c_ndigit_rank_max)-1;

    /**
     * @param ncre
     *  number of creation operators in string
     * @param nann
     *  number of annihilation operators in string
     * @return
     *  max size_t if too many ops in either string
     *  else: a compact, unique integer representing the pair of string lengths
     */
    static constexpr size_t encode_rank_label(size_t ncre, size_t nann) {
        return ((nann>>c_ndigit_rank_max) || (ncre>>c_ndigit_rank_max)) ? ~0ul : (ncre<<c_ndigit_rank_max) | nann;
    }

    static void decode_rank_label(size_t rank_label, size_t& ncre, size_t& nann) {
        ncre = rank_label>>c_ndigit_rank_max;
        if (ncre>>c_ndigit_rank_max){
            // label is out of range
            ncre = ~0ul; nann = ~0ul;
        };
        nann = rank_label&c_ndigit_rank_mask;
    }

    /**
     * indexes all implemented excitation levels and H term ranks
     */
    enum RankIndex {
        rank_00, rank_11, rank_22, rank_01, rank_12, rank_02, rank_10, rank_21, rank_20, rank_count
    };

    /**
     * indexes all contributions to H element for excitation level given by first two digits
     * from H term ranks given by last two digits. See class documentation for more details
     */
    enum ContribCase {
        contrib_0000, contrib_0011, contrib_0022,
        contrib_1111, contrib_1122,
        contrib_2222,
        contrib_0101, contrib_0112,
        contrib_1212,
        contrib_1202,
        contrib_1010, contrib_1021,
        contrib_2121,
        contrib_2020,
        contrib_count
    };

    struct Rank {
        const size_t m_ncre;
        const size_t m_nann;
    };
    struct Contrib {
        const Rank m_exlvl;
        const Rank m_term;
    };

    static constexpr std::array<Contrib, contrib_count> c_contrib_inds = {
            {
                    {{0, 0}, {0, 0}},
                    {{0, 0}, {1, 1}},
                    {{0, 0}, {2, 2}},
                    {{1, 1}, {1, 1}},
                    {{1, 1}, {2, 2}},
                    {{2, 2}, {2, 2}},
                    {{0, 1}, {0, 1}},
                    {{0, 1}, {1, 2}},
                    {{1, 2}, {1, 2}},
                    {{1, 2}, {0, 2}},
                    {{1, 0}, {1, 0}},
                    {{1, 0}, {2, 1}},
                    {{2, 1}, {2, 1}},
                    {{2, 0}, {2, 0}},
            }
    };


    static ContribCase get_coupling_contrib_case(const defs::inds &inds){
        if (inds[2] != ~0ul) {
            // two-body integral
            auto nunique_ind_pairs =
                    (inds[0] != inds[2] && inds[0] != inds[3]) + (inds[1] != inds[2] && inds[1] != inds[3]);
            switch (nunique_ind_pairs) {
                case 0:
                    return contrib_0022;
                case 1:
                    return contrib_1122;
                default:
                    return contrib_2222;
            }
        } else if (inds[0] != ~0ul) {
            // one-body integral
            return inds[0] == inds[1] ? contrib_0011 : contrib_1111;
        } else return contrib_0000;
    }

};


#endif //M7_HAMILTONIANDATA_H











