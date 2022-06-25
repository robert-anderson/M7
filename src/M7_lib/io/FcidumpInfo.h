//
// Created by rja on 25/06/22.
//

#ifndef M7_FCIDUMPINFO_H
#define M7_FCIDUMPINFO_H

#include "FortranNamelistReader.h"
#include "HDF5Wrapper.h"


struct FcidumpInfo {
    const std::string m_fname;
    const bool m_uhf, m_relativistic, m_spin_resolved;
    const uint_t m_nelec, m_nsite, m_nspinorb, m_norb_distinct;
    const int m_ms2;
    const uintv_t m_orbsym;
    FcidumpInfo(std::string fname, bool uhf, bool relativistic, uint_t nelec, uint_t nsite, int ms2, uintv_t orbsym);
    FcidumpInfo(const FortranNamelistReader& reader);
    FcidumpInfo(const hdf5::FileReader& reader);

    FcidumpInfo(std::string fname);
};

struct EbdumpInfo : FcidumpInfo {
    const uint_t m_nmode;
    EbdumpInfo(const FortranNamelistReader& reader):
            FcidumpInfo(reader), m_nmode(reader.read_int("NMODE")){}
    EbdumpInfo(std::string fname): EbdumpInfo(FortranNamelistReader(fname)){}
};

#endif //M7_FCIDUMPINFO_H
