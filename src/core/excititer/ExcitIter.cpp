//
// Created by rja on 25/08/2021.
//

#include "ExcitIter.h"

void ExcitIter::reset(bool need, bool nonzero) {
    m_work_orbs.clear();
    m_need_helement = need;
    m_nonzero_helement_only = nonzero;
    /*
     * even if we don't actually want the helement to be passed to the body function by value, it still needs to be
     * computed so the zero values can be excluded
     */
    if (m_nonzero_helement_only) m_need_helement = true;
}

ExcitIter::ExcitIter(const Hamiltonian &ham, size_t exsig) :
        m_exsig(exsig), m_ham(ham), m_bd(m_ham.m_bd), m_work_conn(m_bd),
        m_work_dst(m_bd), m_work_orbs(ham.m_frm->m_point_group_map) {}
