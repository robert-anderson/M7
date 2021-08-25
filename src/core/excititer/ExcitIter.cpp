//
// Created by rja on 25/08/2021.
//

#include "ExcitIter.h"

void ExcitIter::set_helement_reqs(bool need, bool nonzero) {
    m_need_helement = need;
    m_nonzero_helement_only = nonzero;
    /*
     * even if we don't actually want the helement to be passed to the body function by value, it still needs to be
     * computed so the zero values can be excluded
     */
    if (m_nonzero_helement_only) m_need_helement = true;
}

ExcitIter::ExcitIter(size_t exsig, const Hamiltonian &ham) :
        m_exsig(exsig), m_ham(ham), m_work_conn(ham.nsite()), m_work_dst(ham.nsite()),
        m_work_orbs(ham.m_frm.m_point_group_map){}
