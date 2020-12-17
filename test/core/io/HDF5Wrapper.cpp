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
        auto grp = f.subgroup("first level");
        auto grp2 = grp.subgroup("second level");
        grp2.save(chk, "a complex double");
    }
    {
        hdf5::File f("test.h5", 0);
        auto grp = f.subgroup("first level");
        auto grp2 = grp.subgroup("second level");
        grp2.load(tmp, "a complex double");
        ASSERT_FLOAT_EQ(tmp.real(), chk.real());
        ASSERT_FLOAT_EQ(tmp.imag(), chk.imag());
    }
}


TEST(HDF5Wrapper, Vector){
    hdf5::File f("test.h5", 1);
    auto grp = f.subgroup("container");
    hdf5::VectorWriter<float> vw(grp.m_handle, "my_vector", 10, 1, 10);
    std::vector<float> v = {5, 1, 4, 54, 13524, 3, 2134, 43.12, 123.1243, 456};
    vw.write(v.data(), 1);
}
