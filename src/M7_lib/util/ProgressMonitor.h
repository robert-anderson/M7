//
// Created by Robert J. Anderson on 25/05/2021.
//

#ifndef M7_PROGRESSMONITOR_H
#define M7_PROGRESSMONITOR_H

#include <M7_lib/io/Logging.h>

struct ProgressMonitor {
    const bool m_local;
    const std::string m_name, m_item_name;
    const uint_t m_nexpect, m_pc_resolution, m_period;
    uint_t m_i = 0ul;

    ProgressMonitor(bool local, std::string name, std::string item_name, uint_t nexpect, uint_t pc_resolution = 5);

private:
    void log(const uint_t &pc) const;

public:
    void next();
};


#endif //M7_PROGRESSMONITOR_H
