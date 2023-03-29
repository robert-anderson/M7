//
// Created by Robert J. Anderson on 02/07/2020.
//

#include "Epoch.h"
#include <M7_lib/io/Logging.h>

#include <utility>

Epoch::Epoch(str_t name): m_name(std::move(name)) {
    m_icycle_start.m_local = ~0ul;
    m_icycle_start.m_reduced = ~0ul;
}

bool Epoch::update(uint_t icycle, bool condition) {
    if (*this) return false;
    DEBUG_ASSERT_EQ(icycle_start(), ~0ul, "Epoch should not have already started");
    m_icycle_start.m_local = condition ? icycle : ~0ul;
    m_icycle_start.all_min();
    if (*this) {
        logging::info("Entering \"{}\" epoch on cycle {} ", m_name, icycle);
        return true;
    }
    return false;
}

void Epoch::start(uint_t icycle) {
    update(icycle, true);
}

void Epoch::terminate(uint_t icycle) {
    if (!*this) return;
    ASSERT(icycle_start() != ~0ul)
    logging::info("Terminating \"{}\" epoch on cycle {} ", m_name, icycle);
    m_icycle_start.m_local = ~0ul;
    m_icycle_start.m_reduced = ~0ul;
}

Epoch::operator bool() const {
    return icycle_start() != ~0ul;
}

uint_t Epoch::icycle_start() const {
    return m_icycle_start.m_reduced;
}

bool Epoch::started_last_cycle(uint_t icycle) const {
    return (m_icycle_start.m_reduced!=~0ul && m_icycle_start.m_reduced+1==icycle);
}

bool Epoch::started_this_cycle(uint_t icycle) const {
    return (m_icycle_start.m_reduced!=~0ul && m_icycle_start.m_reduced==icycle);
}

Epochs::Epochs(str_t name, uint_t n, str_t element_identifier) {
    m_epochs.reserve(n);
    for (uint_t i=0ul; i<n; ++i) m_epochs.emplace_back(name+" (" +element_identifier + " " + std::to_string(i) + ")");
}

uint_t Epochs::icycle_start_last() const {
    uint_t highest = 0ul;
    for (const auto& epoch : m_epochs) {
        if (epoch.icycle_start()>highest) highest=epoch.icycle_start();
    }
    return highest;
}

Epochs::operator bool() const {
    return icycle_start_last() != ~0ul;
}