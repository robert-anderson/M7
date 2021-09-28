//
// Created by rja on 15/07/2020.
//

#ifndef M7_FCIDUMPFILEREADER_H
#define M7_FCIDUMPFILEREADER_H

#include "HamiltonianFileReader.h"
#include "src/core/integrals/Integrals_1e.h"
#include "src/core/io/Logging.h"

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

struct FcidumpFileReader : public HamiltonianFileReader {
    const bool m_spin_major;
    const size_t m_nelec;
    const defs::inds m_orbsym;

    bool m_spin_conserving_1e = true;
    bool m_spin_conserving_2e = true;
    size_t m_isymm, m_int_2e_rank;

    using base_t::m_complex_valued;
    FcidumpFileReader(const std::string &fname, bool spin_major);

    bool spin_conserving() const;
    void inds_to_orbs(defs::inds& inds);

    void set_symm_and_rank(const std::string &filename);

    size_t ranksig(const defs::inds &inds) const override;

    size_t exsig(const defs::inds &inds, const size_t& ranksig) const override;

    bool inds_in_range(const defs::inds &inds) const override;
};

#endif //M7_FCIDUMPFILEREADER_H
