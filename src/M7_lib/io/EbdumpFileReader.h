//
// Created by Robert J. Anderson on 22/08/2021.
//

#ifndef M7_EBDUMPFILEREADER_H
#define M7_EBDUMPFILEREADER_H

#include "FcidumpFileReader.h"

struct EbdumpInfo : FcidumpInfo {
    const size_t m_nmode;
    EbdumpInfo(const FortranNamelistReader& reader):
        FcidumpInfo(reader), m_nmode(reader.read_int("NMODE")){}
    EbdumpInfo(std::string fname): EbdumpInfo(FortranNamelistReader(fname)){}
};

struct EbdumpFileReader : HamiltonianFileReader {
    const EbdumpInfo m_info;
    const size_t m_norb_distinct;
    const bool m_spin_major;

    EbdumpFileReader(const std::string &fname, bool spin_major=false);

    size_t ranksig(const defs::ivec_t &inds) const override;

    size_t exsig(const defs::ivec_t &inds, size_t /*ranksig*/) const override;

    bool inds_in_range(const defs::ivec_t &inds) const override;

    void convert_inds(defs::ivec_t &inds);

    bool next(defs::ivec_t &inds, defs::ham_t &v);

};


#endif //M7_EBDUMPFILEREADER_H
