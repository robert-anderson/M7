//
// Created by Robert J. Anderson on 15/07/2020.
//

#ifndef M7_FCIDUMPTEXTFILEREADER_H
#define M7_FCIDUMPTEXTFILEREADER_H

#include <M7_lib/defs.h>
#include <M7_lib/io/Logging.h>

#include "HamTextFileReader.h"
#include "FortranNamelistReader.h"
#include "FcidumpInfo.h"

struct FcidumpTextFileReader : public HamTextFileReader {
    /**
     * TODO: rework this comment to refer to Molpro-style blocks UR style
     * spin-resolved FCIDUMPs index in spinorbs, which may not or may not be spin-major, depending on the program they
     * were generated for. E.g. NECI assumes spin-minor ordering, so if the FCIDUMP supplied was intended for use with
     * NECI, m_ur_style should be false.
     */
    const FcidumpInfo m_info;

    bool m_spin_conserving_1e = true;
    bool m_spin_conserving_2e = true;
private:
    /**
     * each time line of the form "0.0  0 0 0 0" is encountered, increment this counter
     */
    uint_t m_nnull_lines = 0ul;

public:
    FcidumpTextFileReader(const FcidumpInfo& info);

    ~FcidumpTextFileReader() override;

    bool spin_conserving() const;

    /**
     * the index convention employed in the file may differ from that of the coefficient storage, so convert them if
     * necessary
     * @param inds
     *  indices to convert to the storage convention
     */
    void convert_inds(uintv_t &inds);

    bool next(uintv_t &inds, ham_t &v);

    void reset(uint_t iline=0ul);

    OpSig ranksig(const uintv_t &inds) const override;

    OpSig exsig(const uintv_t &inds, OpSig ranksig) const override;

    bool inds_in_range(const uintv_t &inds) const override;

};

#endif //M7_FCIDUMPTEXTFILEREADER_H
