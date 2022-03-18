//
// Created by rja on 02/07/2020.
//

#include "Epoch.h"
#include "io/Logging.h"

Epoch::Epoch(std::string name): m_name(std::move(name)) {
    m_icycle_start.m_local = ~0ul;
    m_icycle_start.m_reduced = ~0ul;
}

bool Epoch::update(size_t icycle, bool condition) {
    if (*this) return false;
    ASSERT(icycle_start() == ~0ul)
    m_icycle_start.m_local = condition?icycle:~0ul;
    m_icycle_start.all_min();
    if (*this) {
        log::info("Entering \"{}\" epoch on cycle {} ", m_name, icycle);
        return true;
    }
    return false;
}

void Epoch::terminate(size_t icycle) {
    if (!*this) return;
    ASSERT(icycle_start() != ~0ul)
    log::info("Terminating \"{}\" epoch on cycle {} ", m_name, icycle);
    m_icycle_start.m_local = ~0ul;
    m_icycle_start.m_reduced = ~0ul;
}

Epoch::operator bool() const {
    return icycle_start() != ~0ul;
}

const size_t &Epoch::icycle_start() const {
    return m_icycle_start.m_reduced;
}

bool Epoch::started_last_cycle(size_t icycle) const {
    return (m_icycle_start.m_reduced!=~0ul && m_icycle_start.m_reduced+1==icycle);
}

bool Epoch::started_this_cycle(size_t icycle) const {
    return (m_icycle_start.m_reduced!=~0ul && m_icycle_start.m_reduced==icycle);
}

size_t Epochs::icycle_start_last() const {
    size_t highest = 0ul;
    for (const auto& epoch : m_epochs) {
        if (epoch.icycle_start()>highest) highest=epoch.icycle_start();
    }
    return highest;
}

Epochs::operator bool() const {
    return icycle_start_last() != ~0ul;
}
