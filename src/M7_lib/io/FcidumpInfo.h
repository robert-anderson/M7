//
// Created by rja on 25/06/22.
//

#ifndef M7_FCIDUMPINFO_H
#define M7_FCIDUMPINFO_H

#include <utility>

#include "M7_lib/hdf5/File.h"
#include "FortranNamelistReader.h"

/**
 * this only involves reading the header of the FCIDUMP file, so it can be done on all ranks - no need for comms
 */
struct FcidumpInfo {
    enum Implementation {
        CSV, MolcasHDF5
    };
    const str_t m_fname;
    const Implementation m_impl;
    const bool m_spin_major;
    /*
     * metadata fields
     */
    bool m_uhf, m_relativistic, m_spin_resolved;
    uint_t m_nelec, m_nsite, m_nspinorb, m_norb_distinct;
    int m_ms2;
    uintv_t m_orbsym;

private:
    void set_metadata(bool uhf, bool relativistic, uint_t nelec, uint_t nsite, int ms2, const uintv_t& orbsym);

public:

    FcidumpInfo(str_t fname, bool spin_major=false);

};

struct EbdumpInfo : FcidumpInfo {
    uint_t m_nmode;
    EbdumpInfo(str_t fname, bool spin_major=false);

};

#endif //M7_FCIDUMPINFO_H
