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

class FcidumpFileReader : public HamiltonianFileReader {
    /**
     * spin-resolved FCIDUMPs index in spinorbs, which may not or may not be spin-major, depending on the program they
     * were generated for. E.g. NECI uses spin-minor ordering throughout, so if the FCIDUMP supplied was intended for
     * use with NECI, spin_major should be passed in as false.
     */
    // spin major and spin restricted (non-resolved) cases
    static void decrement_inds(defs::inds& inds);
    // spin minor case
    static void decrement_inds_and_transpose(defs::inds& inds, const size_t& nspatorb);
    std::function<void(defs::inds& inds)> m_inds_to_orbs;

public:
    const bool m_spin_major;
    const size_t m_nelec;
    const defs::inds m_orbsym;

    bool m_spin_conserving_1e = true;
    bool m_spin_conserving_2e = true;
    size_t m_isymm, m_int_2e_rank;

    using base_t::m_complex_valued;
    FcidumpFileReader(const std::string &fname, bool spin_major);

    bool next(defs::inds &inds, defs::ham_t &v) const override;

    static size_t nind(const defs::inds& inds);

    bool spin_conserving() const;
    void inds_to_orbs(defs::inds& inds);

    void set_symm_and_rank(const std::string &filename);
};

#endif //M7_FCIDUMPFILEREADER_H
