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
        hdf5::FileReader fr("tmp.h5");
        ASSERT_EQ(fr.load_attr<uint_t>("example_attr"), 5ul);
    }
}

TEST(Attributes, ComplexScalar) {
    std::complex<double> z = {1.23, -4.56};
    {
        hdf5::FileWriter fw("tmp.h5");
        fw.save_attr("example_attr", z);
    }
    {
        hdf5::FileReader fr("tmp.h5");
        ASSERT_EQ(fr.load_attr<std::complex<double>>("example_attr"), z);
    }
}

TEST(Attributes, RealVector) {
    const uint_t nitem = 23;
    const auto vec = hash::in_range(123, nitem, 0, 100);
    {
        hdf5::FileWriter fw("tmp.h5");
        fw.save_attr("example_attr", vec);
    }
    {
        hdf5::FileReader fr("tmp.h5");
        ASSERT_EQ(fr.load_attr<v_t<hash::digest_t>>("example_attr"), vec);
    }
}

TEST(Attributes, ComplexVector) {
    const uint_t nitem = 23;
    const auto vec = arith::zip(hash::in_range(123, nitem, 0, 100), hash::in_range(567, nitem, 0, 100));
    {
        hdf5::FileWriter fw("tmp.h5");
        fw.save_attr("example_attr", vec);
    }
    {
        hdf5::FileReader fr("tmp.h5");
        ASSERT_EQ(fr.load_attr<v_t<std::complex<hash::digest_t>>>("example_attr"), vec);
    }
}