//
// Created by Robert John Anderson on 2020-01-17.
//


#include <gtest/gtest.h>
#include "src/core/integrals/Integrals_1e.h"

TEST(Integrals_1e, TwoFoldCheckReal) {
    /*
     * check that all integrals are properly stored by retrieving all
     * values specified in the FCIDUMP and validating against stored values
     */
    typedef double T;
    std::string fname = defs::assets_root+"/RHF_Cr2_12o12e/FCIDUMP";
    Integrals_1e<T, 2> ints(fname);
    FcidumpFileIterator<T> file_iterator(fname);
    defs::inds inds(4);
    T value;
    while (file_iterator.next(inds, value)) {
        if (ints.valid_inds(inds)) {
            ASSERT_TRUE(consts::floats_equal(value, ints.get(inds[0], inds[1])));
        }
    }
}

TEST(Integrals_1e, TwoFoldCheckComplex) {
    /*
     * check that all integrals are properly stored by retrieving all
     * values specified in the FCIDUMP and validating against stored values
     */
    typedef std::complex<double> T;
    std::string fname = defs::assets_root+"/DHF_Be_STO-3G/FCIDUMP";
    Integrals_1e<T, 2> ints(fname);
    FcidumpFileIterator<T> file_iterator(fname);
    defs::inds inds(4);
    T value;
    while (file_iterator.next(inds, value)) {
        if (ints.valid_inds(inds)) {
            ASSERT_TRUE(consts::floats_equal(value, ints.get(inds[0], inds[1])));
        }
    }
}
