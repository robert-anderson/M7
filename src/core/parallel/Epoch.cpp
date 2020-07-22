//
// Created by rja on 02/07/2020.
//

#include "Epoch.h"


Epoch::Epoch(std::string name) :m_name(std::move(name)) {
    m_icycle_start = ~0ul;
    m_icycle_start.mpi_bcast(0);
}

bool Epoch::update(size_t icycle, bool condition) {
    if (*this) return false;
    ASSERT(start() == ~0ul)
    if (condition) {
        m_icycle_start = icycle;
        m_icycle_start.mpi_min();
        ASSERT(start() != ~0ul)
        std::cout << "Entering \"" << m_name << "\" epoch on MC cycle " << icycle << std::endl;
    }
    return *this;
}

Epoch::operator bool() const {
    return start() != ~0ul; }

const size_t &Epoch::start() const {
    return m_icycle_start.reduced();
}
