//
// Created by Robert J. Anderson on 19/08/2021.
//

#ifndef M7_HAMTEXTFILEREADER_H
#define M7_HAMTEXTFILEREADER_H

#include "CsvFileReader.h"
#include "FortranNamelistReader.h"

/**
 * base class for all Hamiltonian-defining files with a Fortran namelist header, and all values defined at the beginning
 * of the line with indices following
 */
class HamTextFileReader : public NumericCsvFileReader {

    strv_t m_work_tokens;

public:
    const bool m_complex_valued;

    static uint_t iline_fn(const str_t& fname);

    HamTextFileReader(const str_t &fname, uint_t nind);

    bool next(uintv_t &inds, ham_t &v);

    static uint_t nset_ind(const uintv_t &inds);

    virtual uint_t ranksig(const uintv_t &inds) const = 0;

    virtual uint_t exsig(const uintv_t &inds, uint_t ranksig) const = 0;

    uint_t exsig(const uintv_t &inds) const;

    virtual bool inds_in_range(const uintv_t& inds) const = 0;

private:
    static void decrement_inds(uintv_t& inds);
};


#endif //M7_HAMTEXTFILEREADER_H
