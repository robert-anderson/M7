//
// Created by Robert John Anderson on 2020-02-09.
//

#include <src/fermion/Determinant.h>
#include "Specification.h"
#include "BitfieldNew.h"

Specification::Specification(const std::array<size_t, ntype> &numeric_lengths,
                             const std::vector<size_t> &bitfield_lengths) :
        m_numeric_lengths(numeric_lengths), m_bitfield_lengths(bitfield_lengths) {
}

Specification::Specification(const std::vector<size_t> &bitfield_lengths) :
        Specification({}, bitfield_lengths) {}

template<>
size_t Specification::add<BitfieldNew>(size_t n) {
    m_bitfield_lengths.push_back(n);
    compile();
    return m_bitfield_lengths.size()-1;
}

template<>
size_t Specification::add<Determinant>(size_t n) {
    /*
     * views on this det will be accessed using the bitfield number of the
     * alpha channel bitfield returned from this function
     */
    add<BitfieldNew>(n); // alpha spin (or Kramers +) channel
    add<BitfieldNew>(n); // beta spin (or Kramers -) channel
    compile();
    return m_bitfield_lengths.size()-2;
}

void Specification::compile() {
    m_numeric_datawords_used = numeric_datawords_used(m_numeric_lengths);
    m_numeric_offsets = numeric_offsets(m_numeric_lengths, m_numeric_datawords_used);
    m_total_numeric_datawords_used = m_numeric_offsets.back() + m_numeric_datawords_used.back();
    m_bitfield_datawords_used = bitfield_datawords_used(m_bitfield_lengths);
    m_bitfield_offsets = bitfield_offsets(m_bitfield_lengths, m_bitfield_datawords_used);
    m_total_bitfield_datawords_used = m_bitfield_lengths.empty() ? 0 :
                                      m_bitfield_offsets.back() + m_bitfield_datawords_used.back();
    m_total_datawords_used = m_total_numeric_datawords_used + m_total_bitfield_datawords_used;
}

std::array<size_t, ntype> Specification::numeric_datawords_used(const std::array<size_t, ntype> &lengths) {
    std::array<size_t, ntype> out{};
    for (size_t i=0ul; i < ntype; ++i) {
        out[i] = (lengths[i] * sizes[i]) / sizeof(defs::data_t) +
                 (((lengths[i] * sizes[i]) % sizeof(defs::data_t)) != 0);
    }
    return out;
}

std::array<size_t, ntype> Specification::numeric_offsets(const std::array<size_t, ntype> &lengths,
                                                         const std::array<size_t, ntype> &datawords_used) {
    std::array<size_t, ntype> out{};
    out[0] = 0;
    for (size_t i = 1ul; i < ntype; ++i) {
        out[i] = out[i - 1] + datawords_used[i - 1];
    }
    return out;
}

std::vector<size_t> Specification::bitfield_datawords_used(const std::vector<size_t> &lengths) {
    std::vector<size_t> out(lengths.size(), 0ul);
    for (size_t i=0ul; i < lengths.size(); ++i) {
        out[i] = lengths[i] / (sizeof(defs::data_t) * 8) + ((lengths[i] % (sizeof(defs::data_t) * 8)) != 0);
    }
    return out;
}

std::vector<size_t>
Specification::bitfield_offsets(const std::vector<size_t> &lengths, const std::vector<size_t> &datawords_used) {
    std::vector<size_t> out(lengths.size(), 0ul);
    for (size_t i=1ul; i < lengths.size(); ++i) {
        out[i] = out[i - 1] + datawords_used[i - 1];
    }
    return out;
}

bool Specification::operator==(const Specification &rhs) const {
    return m_numeric_lengths == rhs.m_numeric_lengths &&
           m_numeric_datawords_used == rhs.m_numeric_datawords_used &&
           m_numeric_offsets == rhs.m_numeric_offsets &&
           m_bitfield_lengths == rhs.m_bitfield_lengths &&
           m_bitfield_datawords_used == rhs.m_bitfield_datawords_used &&
           m_bitfield_offsets == rhs.m_bitfield_offsets &&
           m_total_numeric_datawords_used == rhs.m_total_numeric_datawords_used &&
           m_total_bitfield_datawords_used == rhs.m_total_bitfield_datawords_used &&
           m_total_datawords_used == rhs.m_total_datawords_used;
}

bool Specification::operator!=(const Specification &rhs) const {
    return !(rhs == *this);
}
