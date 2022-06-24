//
// Created by Robert J. Anderson on 19/08/2021.
//

#ifndef M7_HAMILTONIANFILEREADER_H
#define M7_HAMILTONIANFILEREADER_H

#include "CsvFileReader.h"
#include "FortranNamelistReader.h"

/**
 * base class for all Hamiltonian-defining files with a Fortran namelist header, and all values defined at the beginning
 * of the line with indices following
 */
class HamiltonianFileReader : public NumericCsvFileReader {

    std::vector<std::string> m_work_tokens;

public:
    const bool m_complex_valued;

    static uint_t iline_fn(const std::string& fname);

    HamiltonianFileReader(const std::string &fname, uint_t nind);

    bool next(defs::uintv_t &inds, defs::ham_t &v);

    static uint_t nset_ind(const defs::uintv_t &inds);

    virtual uint_t ranksig(const defs::uintv_t &inds) const = 0;

    virtual uint_t exsig(const defs::uintv_t &inds, uint_t ranksig) const = 0;

    uint_t exsig(const defs::uintv_t &inds) const;

    virtual bool inds_in_range(const defs::uintv_t& inds) const = 0;

private:
    static void decrement_inds(defs::uintv_t& inds);
};


#endif //M7_HAMILTONIANFILEREADER_H
