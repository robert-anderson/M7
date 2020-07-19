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


template<typename T>
class FcidumpFileReader : public SparseArrayFileReader<T> {
    /**
     * spin-resolved FCIDUMPs index in spinorbs, which may not may not be spin-major,
     * depending on the program they were generated for. E.g. NECI uses spatial-major
     * ordering throughout, so if the FCIDUMP supplied was intended for use with NECI,
     * spin_major should be passed in as false.
     */
    const bool m_spin_major;
    const size_t m_norb;
    const size_t m_isymm;
    const size_t m_nelec;
    const defs::inds m_orbsym;
    const bool m_spin_resolved;
    const size_t m_nspatorb;
    std::function<void(defs::inds& inds)> m_inds_to_orbs;
    bool m_spin_conserving = true;

    // spin major and spin restricted (non-resolved) cases
    static void decrement_inds(defs::inds& inds){
        for (auto& i:inds) i = ((i==0 || i==~0ul) ? ~0ul : i-1);
    }
    // spin minor case
    static void decrement_inds_and_transpose(defs::inds& inds, const size_t& nspatorb){
        for (auto& i:inds) i = ((i==0 || i==~0ul) ? ~0ul : ((i-1)/2 + ((i&1ul)?0:nspatorb)));
    }

public:
    FcidumpFileReader(const std::string &fname, bool spin_major) : SparseArrayFileReader<T>(fname, 4),
    m_spin_major(spin_major),
    m_norb(read_header_int(fname, "NORB")),
    m_isymm(m_norb>3?isymm(fname):1),
    m_nelec(read_header_int(fname, "NELEC")),
    m_orbsym(read_header_array(fname, "ORBSYM")),
    m_spin_resolved(read_header_bool(fname, "UHF") || read_header_bool(fname, "TREL")),
    m_nspatorb(m_spin_resolved?m_norb/2:m_norb)
    {
        if (m_spin_resolved&!m_spin_major){
            m_inds_to_orbs = [this](defs::inds& inds){
                decrement_inds_and_transpose(inds, m_nspatorb);
            };
        }
        else {
            m_inds_to_orbs = [this](defs::inds& inds){
                decrement_inds(inds);
            };
        }

        if (m_spin_resolved) {
            defs::inds inds(4);
            T v;
            while (next(inds, v)) {
                if (((inds[0]<m_nspatorb)!=(inds[1]<m_nspatorb)) ||
                    ((inds[2]<m_nspatorb)!=(inds[3]<m_nspatorb))){
                    // spin non-conserving example found
                    m_spin_conserving = false;
                    break;
                }
            }
            SparseArrayFileReader<T>::reset(); // go back to beginning of entries
        }
        if (m_spin_conserving) {
            logger::write("FCIDUMP file conserves spin");
        }else{
            logger::write("FCIDUMP file does not conserve spin");
        }
    }

    bool next(defs::inds &inds, T &v) const override {
        auto result = SparseArrayFileReader<T>::next(inds, v);
        m_inds_to_orbs(inds);
        // validate elements
        ASSERT(!result || std::all_of(inds.begin(), inds.end(), [this](const size_t& i){return (i==~0ul)||(i<m_norb);}))
        return result;
    }

    const size_t& norb()const{
        return m_norb;
    }
    const size_t& nelec()const{
        return m_nelec;
    }
    const size_t& nspatorb()const{
        return m_nspatorb;
    }
    const bool& spin_resolved()const{
        return m_spin_resolved;
    }
    const bool& spin_conserving()const{
        return m_spin_conserving;
    }
    void inds_to_orbs(defs::inds& inds){
        m_inds_to_orbs(inds);
    }

    static size_t read_header_int(const std::string &fname, const std::string &label, size_t default_ = 0) {
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
            std::regex_search(line, match, header_terminator_regex);
            if (match.size()) break;
        }
        return default_;
    }

    static defs::inds read_header_array(const std::string &fname, const std::string &label) {
        return defs::inds{0, 1};
    }

    static size_t read_header_bool(const std::string &fname, const std::string &label, size_t default_ = false) {
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
            std::regex_search(line, match, header_terminator_regex);
            if (match.size()) break;
        }
        return default_;
    }


    static size_t isymm(const std::string &filename) {
        std::cout << "Determining permutational symmetry of integral file entries" << std::endl;
        SparseArrayFileReader<T> reader(filename, 4);
        defs::inds inds(4);
        T value;
        // this will eventually hold all orderings of the first example of an
        // integral with 4 distinct indices
        std::array<defs::inds, 8> inds_distinct{};
        size_t isymm = 0ul;
        while (reader.next(inds, value)) {
            if (std::all_of(inds.begin(), inds.end(), [](size_t i){return i>0;})) {
                // we have a two body integral
                if (!isymm) {
                    inds_distinct[0].assign(inds.begin(), inds.end());
                    std::sort(inds.begin(), inds.end());
                    // still looking for an example of four distinct indices
                    if (std::adjacent_find(inds.begin(), inds.end()) == inds.end()) {
                        for (size_t i = 1ul; i < orderings.size(); ++i) {
                            for (size_t j = 0ul; j < 4; ++j)
                                inds_distinct[i].push_back(inds_distinct[0][orderings[i][j]]);
                        }
                        isymm++;
                    }
                } else {
                    for (auto tmp : inds_distinct) {
                        if (std::equal(inds.begin(), inds.end(), tmp.begin())){
                            isymm++;
                            break;
                        }
                    }
                }
            }
        }
        ASSERT(isymm);
        // 8->1; 4->2; 2->4; 1->8
        return 8 / isymm;
    }
};

#endif //M7_FCIDUMPFILEREADER_H
