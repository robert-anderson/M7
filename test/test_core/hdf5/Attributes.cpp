//
// Created by rja on 20/12/22.
//

#include "gtest/gtest.h"
#include "M7_lib/util/Hash.h"
#include "M7_lib/hdf5/File.h"

TEST(Attributes, RealScalar) {
    {
        hdf5::FileWriter fw("tmp.h5");
        fw.save_attr("example_attr", 5ul);
    }
    {
//        hdf5::FileWriter fw("tmp.h5");
//        fw.save_attr("example_attr", 5ul);
    }
}
//TEST(Attributes, RealScalar) {
//    const uint_t nitem = hash::in_range(19 + mpi::irank(), 34, 54);
//    const uint_t max_nitem_per_op = 7;
//    const auto save_vec = hash::in_range(123, nitem, 0, 100);
//    const auto save_vec_all = mpi::all_gatheredv(save_vec);
//    {
//        hdf5::FileWriter fw("tmp.h5");
//        fw.save_dataset("stuff", save_vec, max_nitem_per_op);
//    }
//}