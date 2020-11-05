//
// Created by rja on 05/11/2020.
//

#include "FcidumpFileReader.h"

void FcidumpFileReader::decrement_inds(defs::inds &inds) {
    for (auto& i:inds) i = ((i==0 || i==~0ul) ? ~0ul : i-1);
}

void FcidumpFileReader::decrement_inds_and_transpose(defs::inds &inds, const size_t &nspatorb) {
    for (auto& i:inds) i = ((i==0 || i==~0ul) ? ~0ul : ((i-1)/2 + ((i&1ul)?0:nspatorb)));
}

FcidumpFileReader::FcidumpFileReader(const std::string &fname, bool spin_major) :
        base_t(fname, 4, Tern::False),
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
        defs::ham_t v;
        while (next(inds, v)) {
            if (!consts::float_is_zero(std::abs(v))) {
                if (((inds[0] < m_nspatorb) != (inds[1] < m_nspatorb)) ||
                    ((inds[2] < m_nspatorb) != (inds[3] < m_nspatorb))) {
                    // spin non-conserving example found
                    if (nind(inds)==2) m_spin_conserving_1e = false;
                    else m_spin_conserving_2e = false;
                }
            }
        }
        SparseArrayFileReader<defs::ham_t>::reset(); // go back to beginning of entries
    }
    if (m_spin_conserving_1e) logger::write("FCIDUMP file conserves spin in 1 particle integrals");
    else logger::write("FCIDUMP file does NOT conserve spin in 1 particle integrals");
    if (m_spin_conserving_2e) logger::write("FCIDUMP file conserves spin in 2 particle integrals");
    else logger::write("FCIDUMP file does NOT conserve spin in 2 particle integrals");
}

bool FcidumpFileReader::next(defs::inds &inds, defs::ham_t &v) const {
    auto result = SparseArrayFileReader<defs::ham_t>::next(inds, v);
    m_inds_to_orbs(inds);
    // validate elements
    ASSERT(!result || std::all_of(inds.begin(), inds.end(), [this](const size_t& i){return (i==~0ul)||(i<m_norb);}))
    return result;
}

size_t FcidumpFileReader::nind(const defs::inds &inds) {
    return std::count_if(inds.begin(), inds.end(), [](const size_t& a){return a!=~0ul;});
}

const size_t &FcidumpFileReader::norb() const {
    return m_norb;
}

const size_t &FcidumpFileReader::nelec() const {
    return m_nelec;
}

const size_t &FcidumpFileReader::nspatorb() const {
    return m_nspatorb;
}

const bool &FcidumpFileReader::spin_resolved() const {
    return m_spin_resolved;
}

bool FcidumpFileReader::spin_conserving_1e() const {
    return m_spin_conserving_1e;
}

bool FcidumpFileReader::spin_conserving_2e() const {
    return m_spin_conserving_2e;
}

bool FcidumpFileReader::spin_conserving() const {
    return m_spin_conserving_1e && m_spin_conserving_2e;
}

void FcidumpFileReader::inds_to_orbs(defs::inds &inds) {
    m_inds_to_orbs(inds);
}

size_t FcidumpFileReader::read_header_int(const std::string &fname, const std::string &label, size_t default_) {
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

defs::inds FcidumpFileReader::read_header_array(const std::string &fname, const std::string &label) {
    return defs::inds{0, 1};
}

size_t FcidumpFileReader::read_header_bool(const std::string &fname, const std::string &label, size_t default_) {
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

size_t FcidumpFileReader::isymm(const std::string &filename) {
    std::cout << "Determining permutational symmetry of integral file entries" << std::endl;
    SparseArrayFileReader<defs::ham_t> reader(filename, 4);
    defs::inds inds(4);
    defs::ham_t value;
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
    std::cout << "Permutational symmetry of integral file found to be " << 8/isymm << std::endl;
    // 8->1; 4->2; 2->4; 1->8
    return 8 / isymm;
}
