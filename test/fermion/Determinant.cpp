//
// Created by Robert John Anderson on 2020-02-05.
//

#include <gtest/gtest.h>
#include <src/io/TensorFileIterator.h>

#include "../../src/fermion/Determinant.h"
#include "../../src/enumerators/BitfieldEnumerator.h"
#include "../../src/parallel/MPIWrapper.h"


TEST(Determinant, Phase){

    TensorFileIterator<float> fileIterator(defs::assets_root+"/parity_test/parity_8.txt", 16ul, false);

    defs::inds inds(16);
    float value;

    Determinant bra(4);
    Determinant ket(4);

    while (fileIterator.next(inds, value)){
        bra.zero();
        ket.zero();
        for (auto i{0ul}; i<8; ++i) {
            if (!inds[i]) bra.set(i);
        }
        for (auto i{8ul}; i<16; ++i) {
            if (!inds[i]) ket.set(i-8);
        }
        if (bra.nelec()!=ket.nelec()) continue;
        ASSERT_EQ(bra.phase(ket), ket.phase(bra));
        ASSERT_EQ(bra.phase(ket), value<0);
    }
}