//
// Created by Robert J. Anderson on 25/05/2021.
//

#include "ProgressMonitor.h"

ProgressMonitor::ProgressMonitor(bool local, std::string name, std::string item_name, size_t nexpect,
                                 size_t pc_resolution) :
        m_local(local), m_name(name), m_item_name(item_name), m_nexpect(nexpect),
        m_pc_resolution(pc_resolution), m_period(std::round(nexpect * (pc_resolution / 100.0))) {
    if (m_local)
        log::info_("Starting process \"{}\"...", m_name);
    else
        log::info("Starting process \"{}\"...", m_name);
}

void ProgressMonitor::log(const size_t &pc) const {
    if (m_local)
        log::info_("process \"{}\" is {}% ({}/{} {}) complete",
                   m_name, pc, m_i + 1, m_nexpect, m_item_name);
    else
        log::info("process \"{}\" is {}% ({}/{} {}) complete",
                  m_name, pc, m_i + 1, m_nexpect, m_item_name);
}

void ProgressMonitor::next() {
    auto nperiod = m_i / m_period;
    if (m_i + 1 == m_nexpect) log(100);
    else if (m_i && !(m_i % m_period)) log(nperiod * m_pc_resolution);
    ++m_i;
}
