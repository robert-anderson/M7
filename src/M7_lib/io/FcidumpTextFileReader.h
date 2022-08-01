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

static constexpr std::array<uinta_t<4>, 8> orderings{
        {
                {0, 1, 2, 3},
                {1, 0, 2, 3},
                {0, 1, 3, 2},
                {1, 0, 3, 2},
                {2, 3, 0, 1},
                {2, 3, 1, 0},
                {3, 2, 0, 1},
                {3, 2, 1, 0}
        }
};

struct FcidumpTextFileReader : public HamTextFileReader {
    /**
     * spin-resolved FCIDUMPs index in spinorbs, which may not or may not be spin-major, depending on the program they
     * were generated for. E.g. NECI assumes spin-minor ordering, so if the FCIDUMP supplied was intended for use with
     * NECI, m_spin_major should be false.
     */
    const FcidumpInfo m_info;

    bool m_spin_conserving_1e = true;
    bool m_spin_conserving_2e = true;

    FcidumpTextFileReader(const str_t &fname, bool spin_major);

    bool spin_conserving() const;

    /**
     * the index convention employed in the file may differ from that of the coefficient storage, so convert them if
     * necessary
     * @param inds
     *  indices to convert to the storage convention
     */
    void convert_inds(uintv_t &inds);

    bool next(uintv_t &inds, ham_t &v);

    uint_t ranksig(const uintv_t &inds) const override;

    uint_t exsig(const uintv_t &inds, uint_t ranksig) const override;

    bool inds_in_range(const uintv_t &inds) const override;

};

#endif //M7_FCIDUMPTEXTFILEREADER_H
