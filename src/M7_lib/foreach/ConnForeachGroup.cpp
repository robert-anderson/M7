//
// Created by rja on 12/04/2022.
//

#include "ConnForeachGroup.h"

ConnForeachGroup::ConnForeachGroup(const Hamiltonian &ham) {
    conn_foreach::base_list_t list;
    list = ham.m_frm->make_foreach_iters();
    m_list.merge(list);
    list = ham.m_frmbos->make_foreach_iters();
    m_list.merge(list);
    list = ham.m_bos->make_foreach_iters();
    m_list.merge(list);
    REQUIRE_FALSE(m_list.empty(), "there should be at least one active excitation case");
}
