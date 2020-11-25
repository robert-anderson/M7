//
// Created by rja on 15/07/2020.
//

#ifndef M7_FCIDUMPFILEREADER_H
#define M7_FCIDUMPFILEREADER_H

#include "SparseArrayFileReader.h"
#include <regex>
#include "src/core/integrals/Integrals_1e.h"
#include "src/core/io/Logging.h"

static std::regex header_terminator_regex{R"(\&END)"};

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


class FcidumpFileReader : public SparseArrayFileReader<defs::ham_t> {
    typedef SparseArrayFileReader<defs::ham_t> base_t;
    /**
     * spin-resolved FCIDUMPs index in spinorbs, which may not may not be spin-major,
     * depending on the program they were generated for. E.g. NECI uses spatial-major
     * ordering throughout, so if the FCIDUMP supplied was intended for use with NECI,
     * spin_major should be passed in as false.
     */
    const bool m_spin_major;
    const size_t m_norb;
    const size_t m_nelec;
    const defs::inds m_orbsym;
    const bool m_spin_resolved;
    const size_t m_nspatorb;
    std::function<void(defs::inds& inds)> m_inds_to_orbs;
    bool m_spin_conserving_1e = true;
    bool m_spin_conserving_2e = true;
    size_t m_isymm, m_int_2e_rank;

    // spin major and spin restricted (non-resolved) cases
    static void decrement_inds(defs::inds& inds);
    // spin minor case
    static void decrement_inds_and_transpose(defs::inds& inds, const size_t& nspatorb);

public:
    using base_t::m_complex_valued;
    FcidumpFileReader(const std::string &fname, bool spin_major);

    bool next(defs::inds &inds, defs::ham_t &v) const override;

    static size_t nind(const defs::inds& inds);

    const size_t& norb() const;
    const size_t& nelec() const;
    const size_t& nspatorb() const;
    const bool& spin_resolved() const;
    bool spin_conserving_1e() const;
    bool spin_conserving_2e() const;
    bool spin_conserving() const;
    void inds_to_orbs(defs::inds& inds);

    static size_t read_header_int(const std::string &fname, const std::string &label, size_t default_ = 0);

    static defs::inds read_header_array(const std::string &fname, const std::string &label);

    static size_t read_header_bool(const std::string &fname, const std::string &label, size_t default_ = false);

    void set_symm_and_rank(const std::string &filename);

    const size_t& int_2e_rank() const {return m_int_2e_rank;}
};

#endif //M7_FCIDUMPFILEREADER_H
