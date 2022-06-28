//
// Created by rja on 25/06/22.
//

#ifndef M7_FCIDUMPINFO_H
#define M7_FCIDUMPINFO_H

#include "FortranNamelistReader.h"
#include "M7_lib/hdf5/File.h"


struct FcidumpInfo {
    const std::string m_fname;
    enum Implementation {CSV, MolcasHDF5};
    const Implementation m_impl;
    const bool m_uhf, m_relativistic, m_spin_resolved;
    const uint_t m_nelec, m_nsite, m_nspinorb, m_norb_distinct;
    const int m_ms2;
    const uintv_t m_orbsym;
    FcidumpInfo(std::string fname, Implementation impl,
                bool uhf, bool relativistic, uint_t nelec, uint_t nsite, int ms2, uintv_t orbsym);
    FcidumpInfo(const FortranNamelistReader& reader);
    FcidumpInfo(const hdf5::FileReader& reader);
    /**
     * examine the file at the specified path and decide which ctor to call and return the result
     */
    static FcidumpInfo make(const std::string& fname);
    FcidumpInfo(const std::string& fname): FcidumpInfo(make(fname)){}
};

struct EbdumpInfo : FcidumpInfo {
    const uint_t m_nmode;
    EbdumpInfo(const FortranNamelistReader& reader):
            FcidumpInfo(reader), m_nmode(reader.read_int("NMODE")){}
    EbdumpInfo(std::string fname): EbdumpInfo(FortranNamelistReader(fname)){}
};

#endif //M7_FCIDUMPINFO_H
