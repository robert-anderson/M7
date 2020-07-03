//
// Created by rja on 02/07/2020.
//

#include "Epoch.h"

Epoch::Epoch() {
    m_icycle_start = ~0ul;
    m_icycle_start.mpi_bcast(0);
}

void Epoch::update(size_t icycle, bool condition) {
    if (has_begun()) return;
    if (condition) m_icycle_start = icycle;
    m_icycle_start.mpi_min();
}

const size_t &Epoch::start_cycle() {
    return m_icycle_start.reduced();
}

bool Epoch::has_begun() {
    return start_cycle()!=~0ul;
}
