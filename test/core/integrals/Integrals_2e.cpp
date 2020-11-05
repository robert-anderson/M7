//
// Created by Robert John Anderson on 2020-01-17.
//


#include <gtest/gtest.h>
#include "src/core/integrals/Integrals_2e.h"

TEST(Integrals_2e, FourFoldCheck) {
    /*
     * check that all integrals are properly stored by retrieving all
     * values specified in the FCIDUMP and validating against stored values
     */
    if (!consts::is_complex<defs::ham_t>()) GTEST_SKIP();
    std::string fname = defs::assets_root+"/DHF_Be_STO-3G/FCIDUMP";
    Integrals_2e<defs::ham_t, 4> ints(10, true);
    defs::inds inds(4);
    defs::ham_t value;
    if (mpi::on_node_i_am_root()) {
        FcidumpFileReader file_reader(fname, false);
        while (file_reader.next(inds, value)) {
            if (Integrals_2e<defs::ham_t, 4>::valid_inds(inds)) {
                ints.set(inds, value);
            }
        }
    }
    mpi::barrier_on_node();
    FcidumpFileReader file_reader(fname, false);
    while (file_reader.next(inds, value)) {
        if (Integrals_2e<defs::ham_t, 4>::valid_inds(inds)){
            auto memory_value = ints.element(inds[0], inds[1], inds[2], inds[3]);
            ASSERT_TRUE(consts::floats_equal(value, memory_value));
        }
    }
    //(-0.00233307829512015,-0.00259041061503762)   8  10   6  10  F90 chem spin-minor
    //                                              7   9   5   9  C++ chem spin-minor
    //                                              7   5   9   9  phys spin-minor
    //                                              3   2   4   4  phys spats
    //                                              1   1   1   1  phys spins

    //std::cout << std::endl << ints.get_phys(3, 1, 2, 1, 4, 1, 4, 1) << std::endl;

}
