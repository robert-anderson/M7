//
// Created by rja on 08/05/2021.
//

#include "CiSpaces.h"

#include <utility>

#if 0
ci_gen::Base::Base(size_t nsite, size_t nelec, include_fn_t include_fn) :
    m_nsite(nsite), m_nelec(nelec), m_onv_work(nsite), m_include_fn(std::move(include_fn)){}

ci_gen::NoSym::NoSym(size_t nsite, size_t nelec, const include_fn_t& include_fn) :
        Base(nsite, nelec, include_fn), m_foreach(2*nsite, nelec) {}

void ci_gen::NoSym::operator()(Row &row, fields::Onv<> &onv) {
    ASSERT(onv.belongs_to_row(&row));
    row.m_table->clear();
    auto body = [&](const defs::inds& inds){
        m_onv_work = inds;
        add_if_included(row, onv);
    };
    m_foreach(body);
}

ci_gen::SpinSym::SpinSym(size_t nsite, size_t nelec, int spin, const include_fn_t& include_fn) :
        Base(nsite, nelec, include_fn),
        m_foreach_alpha(nsite, ci_utils::nalpha(nelec, spin)),
        m_foreach_beta(nsite, ci_utils::nbeta(nelec, spin)){}

void ci_gen::SpinSym::operator()(Row &row, fields::FermionOnv &onv) {
    ASSERT(onv.belongs_to_row(&row));
    row.m_table->clear();

    auto alpha_body = [&](const defs::inds& alpha_inds){
        auto beta_body = [&](const defs::inds& beta_inds){
            m_onv_work.zero();
            m_onv_work.set(alpha_inds, beta_inds);
            add_if_included(row, onv);
        };
        m_foreach_beta(beta_body);
    };
    m_foreach_alpha(alpha_body);
}
#endif