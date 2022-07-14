//
// Created by Robert J. Anderson on 25/05/2021.
//

#include "ProgressMonitor.h"

ProgressMonitor::ProgressMonitor(bool local, str_t name, str_t item_name, uint_t nexpect,
                                 uint_t pc_resolution) :
        m_local(local), m_name(name), m_item_name(item_name), m_nexpect(nexpect),
        m_pc_resolution(pc_resolution), m_period(std::round(nexpect * (pc_resolution / 100.0))) {
    if (m_local)
        logging::info_("Starting process \"{}\"...", m_name);
    else
        logging::info("Starting process \"{}\"...", m_name);
}

void ProgressMonitor::log(const uint_t &pc) const {
    if (m_local)
        logging::info_("process \"{}\" is {}% ({}/{} {}) complete",
                   m_name, pc, m_i + 1, m_nexpect, m_item_name);
    else
        logging::info("process \"{}\" is {}% ({}/{} {}) complete",
                  m_name, pc, m_i + 1, m_nexpect, m_item_name);
}

void ProgressMonitor::next() {
    auto nperiod = m_i / m_period;
    if (m_i + 1 == m_nexpect) log(100);
    else if (m_i && !(m_i % m_period)) log(nperiod * m_pc_resolution);
    ++m_i;
}
