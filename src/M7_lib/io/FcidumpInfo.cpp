//
// Created by rja on 25/06/22.
//

#include "FcidumpInfo.h"
#include "M7_lib/basis/BasisData.h"

FcidumpInfo::FcidumpInfo(std::string fname, bool uhf, bool relativistic, uint_t nelec, uint_t nsite, int ms2, uintv_t orbsym):
        m_fname(fname), m_uhf(uhf), m_relativistic(relativistic), m_spin_resolved(m_uhf || m_relativistic),
        m_nelec(nelec), m_nsite(nsite), m_nspinorb(m_spin_resolved ? m_nsite*2 : m_nsite),
        m_norb_distinct(m_spin_resolved ? m_nspinorb : m_nsite),
        m_ms2(ms2), m_orbsym(orbsym.empty() ? uintv_t(m_nsite, 1) : orbsym){
    REQUIRE_EQ(m_orbsym.size(), m_nsite, "invalid ORBSYM specified in FCIDUMP file");
}

FcidumpInfo::FcidumpInfo(const FortranNamelistReader &reader) :
        FcidumpInfo(reader.m_fname, reader.read_bool("UHF"), reader.read_bool("TREL"),
                    reader.read_uint("NELEC"), reader.read_uint("NORB"),
                    reader.read_int("MS2", sys::frm::c_undefined_ms2),
                    integer::dec(reader.read_uints("ORBSYM", {}))){}

FcidumpInfo::FcidumpInfo(const hdf5::FileReader &reader) :
        FcidumpInfo(reader.m_fname,
                    reader.read_attr<int>("UHF", 0),
                    reader.read_attr<int>("TREL", 0),
                    reader.read_attr<uint64_t>("NELEC", 0ul),
                    reader.read_attr<uint64_t>("NORB", 0ul),
                    reader.read_attr<int>("MS2", sys::frm::c_undefined_ms2),
                    integer::dec(reader.read_attr<std::vector<uint64_t>>("ORBSYM", {}))){}

FcidumpInfo::FcidumpInfo(std::string fname) : FcidumpInfo(FortranNamelistReader(fname)){}