//
// Created by anderson on 12/5/21.
//

#ifndef M7_SECTOR_H
#define M7_SECTOR_H

struct Sector {
    // conserved number of electrons in the system
    const size_t m_nelec;
    // maximum occupancy of a single electron mode
    const bool m_elecs;
    // 2 x conserved total z-axis projection of electronic spin (nalpha-nbeta)
    const int m_ms2;
    // conserved number of bosons in the system
    const size_t m_nboson;
    // maximum occupancy of a single boson mode
    const size_t m_nboson_max;
};

#endif //M7_SECTOR_H
