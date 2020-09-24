//
// Created by Robert John Anderson on 2020-03-31.
//

#include "src/core/io/SparseArrayFileReader.h"
#include "gtest/gtest.h"
#include "src/core/basis/Determinant.h"
#include "src/core/basis/Connection.h"

TEST(Determinant, Phase) {

    SparseArrayFileReader<float> file_iterator(defs::assets_root + "/parity_test/parity_8.txt", 16ul);

    defs::inds inds(16);
    float value;

    Determinant bra(4);
    Determinant ket(4);
    AntisymConnection conn(bra);

    while (file_iterator.next(inds, value)) {
        conn.zero();
        bra.zero();
        ket.zero();
        for (auto i{0ul}; i < 8; ++i) {
            if (inds[i]) bra.set(i);
        }
        for (auto i{8ul}; i < 16; ++i) {
            if (inds[i]) ket.set(i - 8);
        }
        if (bra.is_zero() || ket.is_zero()) continue;
        if (bra.nsetbit() != ket.nsetbit()) continue; // TODO: particle number non-conservation
        conn.connect(ket, bra);
        ASSERT_EQ(conn.phase(), value<0);
    }
}

TEST(Determinant, Spin) {
    Determinant det(4);
    ASSERT_EQ(det.spin(), 0);
    det.set(defs::inds{0,1,2,3});
    ASSERT_EQ(det.spin(), 4);
    det.set(defs::inds{4,5,6,7});
    ASSERT_EQ(det.spin(), 0);
    det.zero();
    det.set(defs::inds{4,5,6,7});
    ASSERT_EQ(det.spin(), -4);
    det.zero();
    det.set(defs::inds{0,1,4,5,6,7});
    ASSERT_EQ(det.spin(), -2);
    det.zero();
    det.set(defs::inds{0,1,2,4,5,6,7});
    ASSERT_EQ(det.spin(), -1);
    det.zero();
    det.set(defs::inds{0,1,2,3,4,5,6,7});
    ASSERT_EQ(det.spin(), 0);
    det.zero();
    det.set(defs::inds{0,1,2,3,4,5,6});
    ASSERT_EQ(det.spin(), 1);
}