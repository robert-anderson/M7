//
// Created by rja on 08/05/2021.
//

#include "CiSpaces.h"

ci_gen::Base::Base(size_t nsite, size_t nelec) : m_nsite(nsite), m_nelec(nelec){}

ci_gen::NoSym::NoSym(size_t nsite, size_t nelec) :
        Base(nsite, nelec), m_foreach(2*nsite, nelec) {}

void ci_gen::NoSym::operator()(Row &row, fields::FermionOnv &onv) {
    ASSERT(onv.belongs_to_row(&row));
    row.m_table->clear();
    auto body = [&](const defs::inds& inds){
        row.push_back_jump();
        onv = inds;
        row.m_table->post_insert(row.m_i);
    };
    m_foreach(body);
}

ci_gen::SpinSym::SpinSym(size_t nsite, size_t nelec, int spin) :
        Base(nsite, nelec),
        m_foreach_alpha(nsite, ci_utils::nalpha(nelec, spin)),
        m_foreach_beta(nsite, ci_utils::nbeta(nelec, spin)){}

void ci_gen::SpinSym::operator()(Row &row, fields::FermionOnv &onv) {
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
