//
// Created by Robert John Anderson on 2020-01-17.
//


#include <gtest/gtest.h>
#include "../../src/integrals/Integrals_1e.h"

TEST(Integrals_1e, TwoFoldCheck) {
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
            auto pp = ints.get(inds[0], inds[1]);
            ASSERT_TRUE(consts::floats_equal(value, ints.get(inds[0], inds[1])));
        }
    }
}