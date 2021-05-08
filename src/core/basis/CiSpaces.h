//
// Created by rja on 08/05/2021.
//

#ifndef M7_CISPACES_H
#define M7_CISPACES_H


#include <src/core/enumerator/Enumerator.h>
#include <src/core/field/Fields.h>

namespace ci_gen {

    struct Base {
        const size_t m_nsite, m_nelec;
        Base(size_t nsite, size_t nelec);
    };

    struct NoSym : Base {
        foreach::rtnd::Ordered<> m_foreach;
        NoSym(size_t nsite, size_t nelec);

        void operator()(Row& row, fields::FermionOnv& onv);
    };

    struct SpinSym : Base {
        foreach::rtnd::Ordered<> m_foreach_alpha;
        foreach::rtnd::Ordered<> m_foreach_beta;
        SpinSym(size_t nsite, size_t nelec, int spin);

        void operator()(Row& row, fields::FermionOnv& onv);
    };
};


#endif //M7_CISPACES_H
