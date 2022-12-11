//
// Created by Robert J. Anderson on 22/08/2021.
//

#ifndef M7_EBDUMPFILEREADER_H
#define M7_EBDUMPFILEREADER_H

#include "FcidumpTextFileReader.h"

struct EbdumpFileReader : HamTextFileReader {
    const EbdumpInfo m_info;
    const uint_t m_norb_distinct;

    EbdumpFileReader(const EbdumpInfo& info);

    OpSig ranksig(const uintv_t &inds) const override;

    OpSig exsig(const uintv_t &inds, OpSig /*ranksig*/) const override;

    bool inds_in_range(const uintv_t &inds) const override;

    void convert_inds(uintv_t &inds);

    bool next(uintv_t &inds, ham_t &v);

};


#endif //M7_EBDUMPFILEREADER_H
