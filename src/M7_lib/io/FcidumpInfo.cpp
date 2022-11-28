//
// Created by rja on 25/06/22.
//

#include "FcidumpInfo.h"
#include "M7_lib/basis/BasisData.h"

void FcidumpInfo::set_metadata(bool uhf, bool relativistic, uint_t nelec, uint_t nsite, int ms2, const uintv_t& orbsym) {
    m_uhf = uhf;
    m_relativistic = relativistic;
    m_spin_resolved = m_uhf || m_relativistic;
    m_nelec = nelec;
    m_nsite = nsite;
    m_nspinorb = m_spin_resolved ? m_nsite*2 : m_nsite;
    m_norb_distinct = m_spin_resolved ? m_nspinorb : m_nsite;
    m_ms2 = ms2;
    m_orbsym = orbsym.empty() ? uintv_t(m_nsite, 1) : orbsym;
    REQUIRE_EQ(m_orbsym.size(), m_nsite, "invalid ORBSYM specified in FCIDUMP file");
}

FcidumpInfo::FcidumpInfo(str_t fname, UnrestrictStyle ur_style) :
        m_fname(std::move(fname)), m_impl(hdf5::FileBase::is_hdf5(m_fname) ? MolcasHDF5 : CSV), m_ur_style(ur_style) {
    if (m_impl==CSV) {
        FortranNamelistReader reader(m_fname);
        set_metadata(
            reader.read_bool("UHF") || reader.read_int("IUHF"), reader.read_bool("TREL"),
            reader.read_uint("NELEC"), reader.read_uint("NORB"),
            reader.read_int("MS2", sys::frm::c_undefined_ms2),
            integer::shifted(reader.read_uints("ORBSYM", {}), false));
        if (ur_style==SpinBlocks)
            REQUIRE_FALSE(m_relativistic, "SpinBlocks format is not compatible with spin non-conservation");
    }
    else {
        hdf5::FileReader reader(m_fname);
        set_metadata(
            reader.read_attr<int>("UHF", 0), reader.read_attr<int>("TREL", 0),
            reader.read_attr<int64_t>("NELEC", 0ul),  reader.read_attr<int64_t>("NORB", 0ul),
            reader.read_attr<int>("MS2", sys::frm::c_undefined_ms2),
            integer::shifted(convert::vector<uint_t>(reader.read_attr<v_t<int64_t>>("ORBSYM", {})), false));
    }
}

EbdumpInfo::EbdumpInfo(str_t fname, UnrestrictStyle ur_style) : FcidumpInfo(fname, ur_style) {
    FortranNamelistReader reader(m_fname);
    m_nmode = reader.read_int("NMODE");
}
