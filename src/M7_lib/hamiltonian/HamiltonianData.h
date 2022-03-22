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

        void set_nonzero(size_t exsig);

        bool is_nonzero(size_t exsig) const;

    };

    struct KramersAttributes {
        bool m_conserving_singles = true;
        bool m_conserving_double = true;

        bool conserving() const;
    };
}

#endif //M7_HAMILTONIANDATA_H
