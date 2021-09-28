//
// Created by rja on 19/08/2021.
//

#include "HamiltonianFileReader.h"

HamiltonianFileReader::HamiltonianFileReader(const std::string &fname, size_t nind, bool spin_major) :
        base_t(fname, nind, false),
        m_header_terminator_regex(R"(\&END)"), m_spin_major(spin_major), m_norb(read_header_int(fname, "NORB")),
        m_spin_resolved(read_header_bool(fname, "UHF") || read_header_bool(fname, "TREL")),
        m_nspatorb(m_spin_resolved ? m_norb / 2 : m_norb) {
    if (m_spin_resolved & !m_spin_major) {
        m_inds_to_orbs = [&](defs::inds &inds) {
            decrement_inds_and_transpose(inds, m_nspatorb);
        };
    } else {
        m_inds_to_orbs = [&](defs::inds &inds) {
            decrement_inds(inds);
        };
    }
}

size_t HamiltonianFileReader::read_header_int(const std::string &fname, const std::string &label, size_t default_) {
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

defs::inds HamiltonianFileReader::read_header_array(
        const std::string &fname, const std::string &label, long offset, defs::inds default_) {
    const auto regex = std::regex(label + R"(\s*\=\s*[0-9\,\s]+)");
    FileReader iterator(fname);
    std::string line;
    std::smatch match;
    while (iterator.next(line)) {
        std::regex_search(line, match, regex);
        if (match.size()) {
            // strip away label so any numeric digits in the label are not picked up in error by parse_int_array
            std::string content = match.str();
            std::regex_search(content, match, std::regex(R"(\=.*)"));
            content = match.str();
            return parse_int_array(content, offset);
        }
        std::regex_search(line, match, m_header_terminator_regex);
        if (match.size()) break;
    }
    return default_;
}

size_t HamiltonianFileReader::read_header_bool(const std::string &fname, const std::string &label, size_t default_) {
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

bool HamiltonianFileReader::next(defs::inds &inds, defs::ham_t &v) const {
    auto result = SparseArrayFileReader<defs::ham_t>::next(inds, v);
    if (result) m_inds_to_orbs(inds);
    // validate elements
    DEBUG_ASSERT_TRUE(!result || inds_in_range(inds), "indices in hamiltonian-defining file are not valid");
    return result;
}

size_t HamiltonianFileReader::nset_ind(const defs::inds &inds) {
    return std::count_if(inds.begin(), inds.end(), [](const size_t &a) { return a != ~0ul; });
}

size_t HamiltonianFileReader::exsig(const defs::inds &inds) const {
    return exsig(inds, ranksig(inds));
}

void HamiltonianFileReader::decrement_inds(defs::inds &inds) {
    for (auto &i:inds) i = ((i == 0 || i == ~0ul) ? ~0ul : i - 1);
}

void HamiltonianFileReader::decrement_inds_and_transpose(defs::inds &inds, const size_t &nspatorb) {
    for (auto &i:inds) i = ((i == 0 || i == ~0ul) ? ~0ul : ((i - 1) / 2 + ((i & 1ul) ? 0 : nspatorb)));
}

defs::inds HamiltonianFileReader::parse_int_array(const std::string& str, long offset) {
    defs::inds result;
    std::smatch match;
    auto it = str.cbegin();
    auto end = str.cend();
    while (it!=end){
        std::regex_search(it, end, match, std::regex(R"([0-9]+)"));
        if (!match.size()) break;
        result.push_back(std::atol(match.str().c_str())+offset);
        it+=match.position()+match.size();
    }
    return result;
}