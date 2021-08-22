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

    HamiltonianFileReader(const std::string &fname, size_t nind, bool spin_major):
    base_t(fname, nind, false),
    m_header_terminator_regex(R"(\&END)"), m_spin_major(spin_major), m_norb(read_header_int(fname, "NORB")),
    m_spin_resolved(read_header_bool(fname, "UHF") || read_header_bool(fname, "TREL")),
    m_nspatorb(m_spin_resolved?m_norb/2:m_norb){
        if (m_spin_resolved&!m_spin_major){
            m_inds_to_orbs = [&](defs::inds& inds){
                decrement_inds_and_transpose(inds, m_nspatorb);
            };
        }
        else {
            m_inds_to_orbs = [&](defs::inds& inds){
                decrement_inds(inds);
            };
        }
    }

    size_t read_header_int(const std::string &fname, const std::string &label, size_t default_=0) {
        const auto regex = std::regex(label + R"(\s*\=\s*[0-9]+)");
        FileReader iterator(fname);
        std::string line;
        std::smatch match;
        while (iterator.next(line)) {
            std::regex_search(line, match, regex);
            if (match.size()) {
                std::string match_string(match.str());
                std::regex_search(match_string, match, std::regex(R"([0-9]+)"));
                return std::atol(match.str().c_str());
            }
            std::regex_search(line, match, m_header_terminator_regex);
            if (match.size()) break;
        }
        return default_;
    }

    defs::inds read_header_array(const std::string &fname, const std::string &label) {
        return defs::inds{0, 1};
    }

    size_t read_header_bool(const std::string &fname, const std::string &label, size_t default_=false) {
        std::regex regex;
        if (default_) regex = std::regex(label + R"(\s?\=\s?\.FALSE\.)");
        else regex = std::regex(label + R"(\s?\=\s?\.TRUE\.)");
        FileReader iterator(fname);
        std::string line;
        std::smatch match;
        std::string match_string;
        while (iterator.next(line)) {
            std::regex_search(line, match, regex);
            if (match.size()) return !default_;
            std::regex_search(line, match, m_header_terminator_regex);
            if (match.size()) break;
        }
        return default_;
    }

    bool next(defs::inds &inds, defs::ham_t &v) const {
        auto result = SparseArrayFileReader<defs::ham_t>::next(inds, v);
        if (result) m_inds_to_orbs(inds);
        // validate elements
        ASSERT(!result || std::all_of(inds.begin(), inds.end(), [this](const size_t& i){return (i==~0ul)||(i<m_norb);}))
        return result;
    }

    static size_t nset_ind(const defs::inds &inds) {
        return std::count_if(inds.begin(), inds.end(), [](const size_t& a){return a!=~0ul;});
    }

    virtual size_t ranksig(const defs::inds &inds) const = 0;

    virtual size_t exsig(const defs::inds &inds, const size_t ranksig) const = 0;

    size_t exsig(const defs::inds &inds) const {
        return exsig(inds, ranksig(inds));
    }

private:
    /**
     * spin-resolved FCIDUMPs index in spinorbs, which may not or may not be spin-major, depending on the program they
     * were generated for. E.g. NECI uses spin-minor ordering throughout, so if the FCIDUMP supplied was intended for
     * use with NECI, spin_major should be passed in as false.
     */
    // spin major and spin restricted (non-resolved) cases
    static void decrement_inds(defs::inds& inds){
        for (auto& i:inds) i = ((i==0 || i==~0ul) ? ~0ul : i-1);
    }
    // spin minor case
    static void decrement_inds_and_transpose(defs::inds& inds, const size_t& nspatorb){
        for (auto& i:inds) i = ((i==0 || i==~0ul) ? ~0ul : ((i-1)/2 + ((i&1ul)?0:nspatorb)));
    }
};


#endif //M7_HAMILTONIANFILEREADER_H
