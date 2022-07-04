//
// Created by Robert J. Anderson on 22/08/2021.
//

#ifndef M7_EBDUMPFILEREADER_H
#define M7_EBDUMPFILEREADER_H

#include "FcidumpTextFileReader.h"

struct EbdumpFileReader : HamTextFileReader {
    const EbdumpInfo m_info;
    const uint_t m_norb_distinct;
    const bool m_spin_major;

    EbdumpFileReader(const str_t &fname, bool spin_major=false);

    uint_t ranksig(const uintv_t &inds) const override;

    uint_t exsig(const uintv_t &inds, uint_t /*ranksig*/) const override;

    bool inds_in_range(const uintv_t &inds) const override;

    void convert_inds(uintv_t &inds);

    bool next(uintv_t &inds, ham_t &v);

};


#endif //M7_EBDUMPFILEREADER_H
