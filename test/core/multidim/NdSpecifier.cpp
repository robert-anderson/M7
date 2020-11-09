//
// Created by RJA on 18/09/2020.
//


/*
struct TrivialGetter {
    typedef size_t accessor_t;
    typedef size_t const_accessor_t;
    accessor_t operator()(const size_t& flat) {
        return flat;
    }
};


TEST(NdSpecifier, Indices2D) {
    const size_t nrow = 12;
    const size_t ncol = 7;
    nd_indices::Specifier<2> nd(nrow, ncol);
    size_t iflat = 0ul;
    for (size_t irow = 0ul; irow < nrow; ++irow) {
        for (size_t icol = 0ul; icol < ncol; ++icol) {
            ASSERT_EQ(nd[irow][icol], iflat);
            ++iflat;
        }
    }
}

TEST(NdSpecifier, Indices0D) {
    nd_indices::Specifier<0> nd();
    size_t iflat = 0ul;
    for (size_t ielem = 0ul; ielem < nelem; ++ielem) {
        ASSERT_EQ(nd[ielem], iflat);
        ++iflat;
    }
}

TEST(NdSpecifier, Indices1D) {
    const size_t nelem = 12;
    nd_indices::Specifier<1> nd(nelem);
    size_t iflat = 0ul;
    for (size_t ielem = 0ul; ielem < nelem; ++ielem) {
        ASSERT_EQ(nd[ielem], iflat);
        ++iflat;
    }
}

TEST(ArrayFormat, NdAccessor) {
    nd_array::Specifier<float, 2> array(3, 4);
    array[0][3] = 2;
    ASSERT_EQ(array.select(3), array[0][3]);

    std::cout << array.to_string() << std::endl;
}

TEST(NdSpecifier, Sequence) {
    NdSpecifier<TrivialGetter, 3> af(TrivialGetter(), 4, 2, 3);
    size_t i = 0ul;
    for (size_t i0 = 0ul; i0 < af.shape()[0]; ++i0) {
        for (size_t i1 = 0ul; i1 < af.shape()[1]; ++i1) {
            for (size_t i2 = 0ul; i2 < af.shape()[2]; ++i2) {
                ASSERT_EQ(af[i0][i1][i2].select(), i);
                ++i;
            }
        }
    }
}
 */


/*
TEST(NdSpecifier, OneDimCase){
    NdSpecifier<1> af(12);
    for (size_t i = 0ul; i < af.nelement(); ++i) {
        ASSERT_EQ(af[i], i);
    }
}

TEST(NdSpecifier, Sequence) {
    NdSpecifier<3> af(4, 2, 3);
    size_t i = 0ul;
    for (size_t i0 = 0ul; i0 < af.shape()[0]; ++i0) {
        for (size_t i1 = 0ul; i1 < af.shape()[1]; ++i1) {
            for (size_t i2 = 0ul; i2 < af.shape()[2]; ++i2) {
                ASSERT_EQ(af[i0][i1][i2], i);
                ++i;
            }
        }
    }
}
 */