//
// Created by rja on 13/12/2020.
//

#include "gtest/gtest.h"
#include "src/core/io/HDF5Wrapper.h"

TEST(HDF5Wrapper, Complex){
    std::complex<double> chk = {1.213, 213.346};
    std::complex<double> tmp;
    {
        hdf5::File f("test.h5", 1);
        auto grp = f.group("first level");
        //auto grp2 = grp.subgroup("second level");
        grp.save(chk, "a complex double");
    }
    {
        hdf5::File f("test.h5", 0);
        auto grp = f.group("first level");
        //auto grp2 = grp.subgroup("second level");
        grp.load(tmp, "a complex double");
        ASSERT_FLOAT_EQ(tmp.real(), chk.real());
        ASSERT_FLOAT_EQ(tmp.imag(), chk.imag());
    }
}