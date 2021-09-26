//
// Created by rja on 19/08/2021.
//

#ifndef M7_HAMILTONIANFILEREADER_H
#define M7_HAMILTONIANFILEREADER_H

#include <regex>
#include "SparseArrayFileReader.h"

/**
 * base class for all Hamiltonian-defining files with a FORTRAN namelist header, and all values defined at the beginning
 * of the line with indices following
 */
struct HamiltonianFileReader : public SparseArrayFileReader<defs::ham_t> {
    typedef SparseArrayFileReader<defs::ham_t> base_t;
    std::function<void(defs::inds& inds)> m_inds_to_orbs;

    const std::regex m_header_terminator_regex;
    const bool m_spin_major;
    const size_t m_norb;
    const bool m_spin_resolved;
    const size_t m_nspatorb;

    HamiltonianFileReader(const std::string &fname, size_t nind, bool spin_major);

    size_t read_header_int(const std::string &fname, const std::string &label, size_t default_=0);

    defs::inds read_header_array(const std::string &fname, const std::string &label, long offset=0, defs::inds default_={});

    size_t read_header_bool(const std::string &fname, const std::string &label, size_t default_=false);

    bool next(defs::inds &inds, defs::ham_t &v) const;

    static size_t nset_ind(const defs::inds &inds);

    virtual size_t ranksig(const defs::inds &inds) const = 0;

    virtual size_t exsig(const defs::inds &inds, const size_t& ranksig) const = 0;

    size_t exsig(const defs::inds &inds) const;

private:
    /**
     * spin-resolved FCIDUMPs index in spinorbs, which may not or may not be spin-major, depending on the program they
     * were generated for. E.g. NECI uses spin-minor ordering throughout, so if the FCIDUMP supplied was intended for
     * use with NECI, spin_major should be passed in as false.
     */
    // spin major and spin restricted (non-resolved) cases
    static void decrement_inds(defs::inds& inds);
    // spin minor case
    static void decrement_inds_and_transpose(defs::inds& inds, const size_t& nspatorb);

    static defs::inds parse_int_array(const std::string& str, long offset=0);
};


#endif //M7_HAMILTONIANFILEREADER_H
