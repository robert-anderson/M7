//
// Created by Robert J. Anderson on 05/11/2020.
//

#include "FcidumpTextFileReader.h"
#include "M7_lib/util/Exsig.h"
#include "M7_lib/basis/BasisData.h"

FcidumpTextFileReader::FcidumpTextFileReader(const std::string &fname, bool spin_major) :
        HamTextFileReader(fname, 4), m_info(FortranNamelistReader(fname)), m_spin_major(spin_major) {
    auto nsite = m_info.m_nsite;
    if (m_info.m_spin_resolved) {
        uintv_t inds(4);
        ham_t v;
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

bool FcidumpTextFileReader::spin_conserving() const {
    return m_spin_conserving_1e && m_spin_conserving_2e;
}

void FcidumpTextFileReader::convert_inds(uintv_t &inds) {
    if (!m_info.m_spin_resolved || m_spin_major) return;
    for (auto &ind: inds)
        ind = (ind == ~0ul) ? ~0ul : (ind / 2 + ((ind & 1ul) ? m_info.m_nsite : 0));
}

bool FcidumpTextFileReader::next(uintv_t &inds, ham_t &v) {
    if (!HamTextFileReader::next(inds, v)) return false;
    convert_inds(inds);
    return true;
}

uint_t FcidumpTextFileReader::ranksig(const uintv_t &inds) const {
    auto nset_inds = HamTextFileReader::nset_ind(inds);
    return exsig::encode(nset_inds/2, nset_inds/2, 0, 0);
}

uint_t FcidumpTextFileReader::exsig(const uintv_t &inds, uint_t ranksig) const {
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

bool FcidumpTextFileReader::inds_in_range(const uintv_t &inds) const {
    return std::all_of(inds.cbegin(), inds.cend(), [this](uint_t i){return i<=m_info.m_norb_distinct;});
}
