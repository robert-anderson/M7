//
// Created by rja on 08/05/2021.
//

#include "CiSpaces.h"

#include <utility>

ci_gen::Base::Base(BasisData bd, size_t nelec, include_fn_t include_fn) :
    m_bd(bd), m_nelec(nelec), m_mbf_work(bd), m_include_fn(std::move(include_fn)){}

ci_gen::NoSym::NoSym(BasisData bd, size_t nelec, const include_fn_t& include_fn) :
        Base(bd, nelec, include_fn), m_foreach(2*bd.m_nsite, nelec) {}

ci_gen::SpinSym::SpinSym(BasisData bd, size_t nelec, int spin, const include_fn_t& include_fn) :
        Base(bd, nelec, include_fn),
        m_foreach_alpha(bd.m_nsite, ci_utils::nalpha(nelec, spin)),
        m_foreach_beta(bd.m_nsite, ci_utils::nbeta(nelec, spin)){}

void ci_gen::SpinSym::operator()(Row &row, field::Mbf &mbf) {
    ASSERT(mbf.belongs_to_row(&row));
    row.m_table->clear();

    auto body = [&](){
        m_mbf_work[mbf].zero();
        set_from_inds(m_mbf_work[mbf], m_foreach_alpha.inds(), m_foreach_beta.inds());
        add_if_included(row, mbf);
    };
    m_foreach_alpha(body, m_foreach_beta);
}
