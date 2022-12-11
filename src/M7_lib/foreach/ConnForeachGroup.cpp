//
// Created by Robert J. Anderson on 12/04/2022.
//

#include "ConnForeachGroup.h"

ConnForeachGroup::ConnForeachGroup(const Hamiltonian &ham) {
    conn_foreach::base_list_t list;
    list = ham.m_frm.make_foreach_iters();
    m_list.merge(list);
    list = ham.m_frmbos.make_foreach_iters();
    m_list.merge(list);
    list = ham.m_bos.make_foreach_iters();
    m_list.merge(list);
    REQUIRE_FALSE(m_list.empty(), "there should be at least one active excitation case");
    log();
}

void ConnForeachGroup::log() const {
    strv_t exsigs;
    for (auto& it : m_list) exsigs.emplace_back(it->m_exsig.to_string());
    logging::info("Connection iterator excitation signatures: {}", convert::to_string(exsigs));
}
