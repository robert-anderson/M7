//
// Created by Robert John Anderson on 2020-03-31.
//

#include "src/core/io/SparseArrayFileReader.h"
#include "gtest/gtest.h"
#include "src/core/field/Elements.h"
#include "src/core/basis/FermionOnvConnection.h"

TEST(Determinant, Phase) {

    SparseArrayFileReader<float> file_reader(
            defs::assets_root + "/parity_test/parity_8.txt",
            16ul, false, false);

    defs::inds inds(16);
    float value;

    const size_t nsite = 4;
    elements::FermionOnv bra(nsite);
    elements::FermionOnv ket(nsite);
    AntisymFermionOnvConnection conn(bra);

    while (file_reader.next(inds, value)) {
        conn.zero();
        bra.zero();
        ket.zero();
        for (size_t i = 0ul; i < 8ul; ++i) {
            if (inds[i]) bra.set(i);
        }
        for (size_t i = 8ul; i < 16ul; ++i) {
            if (inds[i]) ket.set(i - 8);
        }
        if (bra.is_zero() || ket.is_zero()) continue;
        if (bra.nsetbit() != ket.nsetbit()) continue; // TODO: particle number non-conservation
        conn.connect(ket, bra);
        ASSERT_EQ(conn.phase(), value<0);
    }
}

TEST(Determinant, Spin) {
    elements::FermionOnv det(4);
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
