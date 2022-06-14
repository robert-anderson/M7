//
// Created by Robert J. Anderson on 15/07/2020.
//

#ifndef M7_FCIDUMPFILEREADER_H
#define M7_FCIDUMPFILEREADER_H

#include <M7_lib/defs.h>
#include <M7_lib/io/Logging.h>

#include "HamiltonianFileReader.h"
#include "FortranNamelistReader.h"

static constexpr std::array<std::array<size_t, 4>, 8> orderings{
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

struct FcidumpInfo {
    const std::string m_fname;
    const bool m_uhf, m_relativistic, m_spin_resolved;
    const size_t m_nelec, m_nsite, m_nspinorb, m_norb_distinct;
    const int m_ms2;
    const defs::inds m_orbsym;
    FcidumpInfo(std::string fname, bool uhf, bool relativistic, size_t nelec, size_t nsite, int ms2, defs::inds orbsym);
    FcidumpInfo(const FortranNamelistReader& reader);

    FcidumpInfo(std::string fname);
};

struct FcidumpFileReader : public HamiltonianFileReader {
    /**
     * spin-resolved FCIDUMPs index in spinorbs, which may not or may not be spin-major, depending on the program they
     * were generated for. E.g. NECI assumes spin-minor ordering, so if the FCIDUMP supplied was intended for use with
     * NECI, m_spin_major should be false.
     */
    const FcidumpInfo m_info;
    const bool m_spin_major;

    bool m_spin_conserving_1e = true;
    bool m_spin_conserving_2e = true;

    FcidumpFileReader(const std::string &fname, bool spin_major);

    bool spin_conserving() const;

    /**
     * the index convention employed in the file may differ from that of the coefficient storage, so convert them if
     * necessary
     * @param inds
     *  indices to convert to the storage convention
     */
    void convert_inds(defs::inds &inds);

    bool next(defs::inds &inds, defs::ham_t &v);

    size_t ranksig(const defs::inds &inds) const override;

    size_t exsig(const defs::inds &inds, const size_t &ranksig) const override;

    bool inds_in_range(const defs::inds &inds) const override;

};

#endif //M7_FCIDUMPFILEREADER_H
