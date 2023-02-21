//
// Created by Robert J. Anderson on 19/08/2021.
//

#include "HamTextFileReader.h"

uint_t HamTextFileReader::iline_fn(const str_t &fname) {
    FileReader file_reader(fname);
    str_t line;
    uint_t iline=0ul;
    while (file_reader.next(line)){
        ++iline;
        if (FortranNamelistReader::contains_terminator(line)) return iline;
    }
    return ~0ul;
}

HamTextFileReader::HamTextFileReader(const str_t &fname, uint_t nind) :
        NumericCsvFileReader(fname, ncolumn(fname, iline_fn)),
        m_complex_valued(m_ncolumn==nind+2) {
    REQUIRE_GT(m_ncolumn, nind, "not enough columns in body of file "+fname);
    REQUIRE_LE(m_ncolumn, nind+2, "too many columns in body of file "+fname);
    if (!dtype::is_complex<ham_t>()){
        REQUIRE_FALSE(m_complex_valued, "can't read complex-valued array into a real container");
    }
}

bool HamTextFileReader::next(uintv_t &inds, ham_t &v) {
    auto result = NumericCsvFileReader::next(m_work_tokens);
    if (!result) return false;
    REQUIRE_EQ(m_work_tokens.size(), m_ncolumn, "invalid line found in file "+m_fname);
    auto begin = m_work_tokens.cbegin();
    auto ind_begin = begin + (m_complex_valued ? 2 : 1);
    parse::checked(begin, ind_begin, v);
    parse::checked(ind_begin, m_work_tokens.cend(), inds);
    REQUIRE_TRUE(inds_in_range(inds), "index OOB");
    // 1-based indexing is assumed
    decrement_inds(inds);
    return true;
}

uint_t HamTextFileReader::nset_ind(const uintv_t &inds) {
    return std::count_if(inds.begin(), inds.end(), [](const uint_t &a) { return a != ~0ul; });
}

OpSig HamTextFileReader::exsig(const uintv_t &inds) const {
    return exsig(inds, ranksig(inds));
}

void HamTextFileReader::decrement_inds(uintv_t &inds) {
    for (auto &i:inds) i = ((i == 0 || i == ~0ul) ? ~0ul : i - 1);
}