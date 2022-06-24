//
// Created by Robert J. Anderson on 05/11/2020.
//

#include "FcidumpFileReader.h"
#include "M7_lib/util/Exsig.h"
#include "M7_lib/basis/BasisData.h"

FcidumpInfo::FcidumpInfo(std::string fname, bool uhf, bool relativistic, uint_t nelec, uint_t nsite, int ms2, defs::uintv_t orbsym):
        m_fname(fname), m_uhf(uhf), m_relativistic(relativistic), m_spin_resolved(m_uhf || m_relativistic),
        m_nelec(nelec), m_nsite(nsite), m_nspinorb(m_spin_resolved ? m_nsite*2 : m_nsite),
        m_norb_distinct(m_spin_resolved ? m_nspinorb : m_nsite),
        m_ms2(ms2), m_orbsym(orbsym.empty() ? defs::uintv_t(m_nsite, 1) : orbsym){
    REQUIRE_EQ(m_orbsym.size(), m_nsite, "invalid ORBSYM specified in FCIDUMP file");
}

FcidumpInfo::FcidumpInfo(const FortranNamelistReader &reader) :
        FcidumpInfo(reader.m_fname, reader.read_bool("UHF"), reader.read_bool("TREL"),
                    reader.read_uint("NELEC"), reader.read_uint("NORB"),
                    reader.read_int("MS2", sys::frm::c_undefined_ms2),
                    integer::dec(reader.read_uints("ORBSYM", {}))){}

FcidumpInfo::FcidumpInfo(std::string fname) : FcidumpInfo(FortranNamelistReader(fname)){}

FcidumpFileReader::FcidumpFileReader(const std::string &fname, bool spin_major) :
        HamiltonianFileReader(fname, 4), m_info(FortranNamelistReader(fname)), m_spin_major(spin_major) {
    auto nsite = m_info.m_nsite;
    if (m_info.m_spin_resolved) {
        defs::uintv_t inds(4);
        defs::ham_t v;
        while (next(inds, v)) {
            if (fptol::numeric_zero(v)) continue;
            if (((inds[0] < nsite) != (inds[1] < nsite)) || ((inds[2] < nsite) != (inds[3] < nsite))) {
                // spin non-conserving example found
                if (nset_ind(inds)==2) m_spin_conserving_1e = false;
                else m_spin_conserving_2e = false;
            }
        }
        FileReader::reset(); // go back to beginning of entries
    }
    if (m_spin_conserving_1e) log::info("FCIDUMP file conserves spin in 1 particle integrals");
    else log::info("FCIDUMP file does NOT conserve spin in 1 particle integrals");
    if (m_spin_conserving_2e) log::info("FCIDUMP file conserves spin in 2 particle integrals");
    else log::info("FCIDUMP file does NOT conserve spin in 2 particle integrals");
}

bool FcidumpFileReader::spin_conserving() const {
    return m_spin_conserving_1e && m_spin_conserving_2e;
}

void FcidumpFileReader::convert_inds(defs::uintv_t &inds) {
    if (!m_info.m_spin_resolved || m_spin_major) return;
    for (auto &ind: inds)
        ind = (ind == ~0ul) ? ~0ul : (ind / 2 + ((ind & 1ul) ? m_info.m_nsite : 0));
}

bool FcidumpFileReader::next(defs::uintv_t &inds, defs::ham_t &v) {
    if (!HamiltonianFileReader::next(inds, v)) return false;
    convert_inds(inds);
    return true;
}

uint_t FcidumpFileReader::ranksig(const defs::uintv_t &inds) const {
    auto nset_inds = HamiltonianFileReader::nset_ind(inds);
    return exsig::encode(nset_inds/2, nset_inds/2, 0, 0);
}

uint_t FcidumpFileReader::exsig(const defs::uintv_t &inds, uint_t ranksig) const {
    switch (ranksig) {
        case 0ul: return 0ul;
        case exsig::ex_single:
            return inds[0]==inds[1] ? 0ul : exsig::ex_single;
            case exsig::ex_double:
            return inds[0]==inds[2] ?
                (inds[1]==inds[3] ? 0ul : exsig::ex_single):
                (inds[1]==inds[3] ? exsig::ex_single : exsig::ex_double);
        default:
            return ~0ul;
    }
}

bool FcidumpFileReader::inds_in_range(const defs::uintv_t &inds) const {
    return std::all_of(inds.cbegin(), inds.cend(), [this](uint_t i){return i<=m_info.m_norb_distinct;});
}
