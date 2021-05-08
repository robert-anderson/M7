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
        Base(size_t nsite, size_t nelec) : m_nsite(nsite), m_nelec(nelec){}
    };

    struct NoSym : Base {
        foreach::rtnd::Ordered<> m_foreach;
        NoSym(size_t nsite, size_t nelec) :
        Base(nsite, nelec), m_foreach(2*nsite, nelec) {}

        void operator()(Row& row, fields::FermionOnv& onv) {
            ASSERT(onv.belongs_to_row(&row));
            row.m_table->clear();
            auto body = [&](const defs::inds& inds){
                row.push_back_jump();
                onv = inds;
                row.m_table->post_insert(row.m_i);
            };
            m_foreach(body);
        }
    };

    struct SpinSym : Base {
        foreach::rtnd::Ordered<> m_foreach_alpha;
        foreach::rtnd::Ordered<> m_foreach_beta;
        SpinSym(size_t nsite, size_t nelec, int spin) :
                Base(nsite, nelec),
                m_foreach_alpha(nsite, ci_utils::nalpha(nelec, spin)),
                m_foreach_beta(nsite, ci_utils::nbeta(nelec, spin)){}

        void operator()(Row& row, fields::FermionOnv& onv) {
            ASSERT(onv.belongs_to_row(&row));
            row.m_table->clear();

            auto alpha_body = [&](const defs::inds& alpha_inds){
                auto beta_body = [&](const defs::inds& beta_inds){
                    row.push_back_jump();
                    onv.set(alpha_inds, beta_inds);
                    row.m_table->post_insert(row.m_i);
                };
                m_foreach_beta(beta_body);
            };
            m_foreach_alpha(alpha_body);
        }
    };
};


#endif //M7_CISPACES_H
