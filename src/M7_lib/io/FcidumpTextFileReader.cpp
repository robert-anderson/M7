//
// Created by Robert J. Anderson on 05/11/2020.
//

#include "FcidumpTextFileReader.h"
#include "M7_lib/basis/BasisData.h"

FcidumpTextFileReader::FcidumpTextFileReader(const FcidumpInfo& info) :
        HamTextFileReader(info.m_fname, 4), m_info(info){
    auto nsite = m_info.m_nsite;
    if (m_info.m_spin_resolved && (m_info.m_ur_style!=FcidumpInfo::SpinBlocks)) {
        uintv_t inds(4);
        ham_t v;
        while (next(inds, v)) {
            if (fptol::near_zero(v)) continue;
            if (((inds[0] < nsite) != (inds[1] < nsite)) || ((inds[2] < nsite) != (inds[3] < nsite))) {
                // spin non-conserving example found
                if (nset_ind(inds)==2) m_spin_conserving_1e = false;
                else m_spin_conserving_2e = false;
            }
        }
        reset(); // go back to beginning of entries
    }
    logging::info("FCIDUMP file {} spin in 1 particle integrals",
                  m_spin_conserving_1e ? "conserves" : "does NOT conserve");
    logging::info("FCIDUMP file {} spin in 2 particle integrals",
                  m_spin_conserving_2e ? "conserves" : "does NOT conserve");
}

FcidumpTextFileReader::~FcidumpTextFileReader() {
    if ((m_nnull_lines > 1) && m_info.m_ur_style!=FcidumpInfo::SpinBlocks) {
        logging::warn(R"(More than one block delimiter line "0.0  0  0  0  0" found, but "blocks" unrestrict style not enabled)");
    }
}

bool FcidumpTextFileReader::spin_conserving() const {
    return m_spin_conserving_1e && m_spin_conserving_2e;
}

void FcidumpTextFileReader::convert_inds(uintv_t &inds) {
    if (!m_info.m_spin_resolved || m_info.m_ur_style==FcidumpInfo::SpinMajor) return;
    if (m_info.m_spin_resolved && m_info.m_ur_style==FcidumpInfo::SpinBlocks) {
        switch (m_nnull_lines) {
            case 0:
                // (uu|uu) block
                break;
            case 1:
                // (dd|dd) block
                inds[0] += m_info.m_nsite;
                inds[1] += m_info.m_nsite;
                inds[2] += m_info.m_nsite;
                inds[3] += m_info.m_nsite;
                break;
            case 2:
                // (uu|dd) block
                inds[2] += m_info.m_nsite;
                inds[3] += m_info.m_nsite;
                break;
            case 3:
                // h_uu block
                break;
            case 4:
                // h_dd block
                inds[0] += m_info.m_nsite;
                inds[1] += m_info.m_nsite;
                break;
            default:
                // core energy
                REQUIRE_TRUE(std::all_of(inds.cbegin(), inds.cend(), [](uint_t i){return i==~0ul;}),
                             "integrals specified beyond last block");
                break;
        }
        return;
    }
    for (auto &ind: inds)
        ind = (ind == ~0ul) ? ~0ul : (ind / 2 + ((ind & 1ul) ? m_info.m_nsite : 0));
}

bool FcidumpTextFileReader::next(uintv_t &inds, ham_t &v) {
    if (!HamTextFileReader::next(inds, v)) return false;
    // check for line of the form "0.0  0 0 0 0"
    if ((v==0.0) && std::all_of(inds.cbegin(), inds.cend(), [](uint_t i){return i==~0ul;})){
        ++m_nnull_lines;
        return next(inds, v);
    }
    convert_inds(inds);
    return true;
}

void FcidumpTextFileReader::reset(uint_t iline) {
    FileReader::reset(iline);
    m_nnull_lines = 0ul;
}

OpSig FcidumpTextFileReader::ranksig(const uintv_t &inds) const {
    auto nset_inds = HamTextFileReader::nset_ind(inds);
    return opsig::frm(nset_inds/2);
}

OpSig FcidumpTextFileReader::exsig(const uintv_t &inds, OpSig ranksig) const {
    switch (ranksig.to_int()) {
        case opsig::c_0000.to_int(): return opsig::c_0000;
        case opsig::c_sing.to_int():
            return (inds[0]==inds[1]) ? opsig::c_0000 : opsig::c_sing;
            case opsig::c_doub.to_int():
            return inds[0]==inds[2] ?
                (inds[1]==inds[3] ? opsig::c_0000 : opsig::c_sing):
                (inds[1]==inds[3] ? opsig::c_sing : opsig::c_doub);
        default:
            return opsig::c_invalid;
    }
}

bool FcidumpTextFileReader::inds_in_range(const uintv_t &inds) const {
    return std::all_of(inds.cbegin(), inds.cend(), [this](uint_t i){return i<=m_info.m_norb_distinct;});
}
