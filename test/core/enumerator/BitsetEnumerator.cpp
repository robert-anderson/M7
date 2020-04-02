//
// Created by Robert John Anderson on 2020-03-30.
//

#include <src/core/fermion/Determinant.h>
#include <src/core/fermion/DecodedDeterminant.h>
#include "gtest/gtest.h"
#include "src/core/enumerator/BitsetEnumerator.h"
#include "src/core/table/Bitset.h"

TEST(BitsetEnumerator, BitsetSet) {
    Bitset bitset(200);
    defs::inds setbits{1, 4, 7, 32, 89, 123, 199};
    bitset.set(setbits);
    BitsetSetEnumerator enumerator(bitset);

    size_t setbit;
    size_t i = ~0ul;
    while (enumerator.next(setbit, i)) {
        ASSERT_EQ(setbits[i], setbit);
    }
}

TEST(BitsetEnumerator, BitsetClr) {
    Bitset bitset(200);
    defs::inds setbits{1, 4, 7, 32, 89, 123, 199};

    bitset.set(setbits);
    BitsetClrEnumerator enumerator(bitset);

    auto setptr = setbits.begin();
    size_t clrbit;
    size_t i = ~0ul;
    while (enumerator.next(clrbit, i)) {
        if (clrbit > *setptr) setptr++;
        if (setptr == setbits.end()) break;
        ASSERT_NE(*setptr, clrbit);
    }
}

TEST(BitsetEnumerator, BitsetAnd) {
    Bitset bitset1(200);
    defs::inds setbits1{1, 4, 7, 32, 89, 123, 199};
    Bitset bitset2(200);
    defs::inds setbits2{1, 4, 5, 7, 89, 91, 123};

    defs::inds andbits{1, 4, 7, 89, 123};

    bitset1.set(setbits1);
    bitset2.set(setbits2);

    BitsetAndEnumerator enumerator(bitset1, bitset2);

    size_t andbit;
    size_t i = ~0ul;
    while (enumerator.next(andbit, i)) {
        ASSERT_EQ(andbits[i], andbit);
    }
}

TEST(BitsetEnumerator, BitsetAndNot) {
    Bitset bitset1(200);
    defs::inds setbits1{1, 4, 7, 32, 89, 124, 199};
    Bitset bitset2(200);
    defs::inds setbits2{1, 4, 5, 7, 89, 91, 123};

    defs::inds andnotbits{32, 124, 199};

    bitset1.set(setbits1);
    bitset2.set(setbits2);

    BitsetAndNotEnumerator enumerator(bitset1, bitset2);

    size_t andnotbit;
    size_t i = ~0ul;
    while (enumerator.next(andnotbit, i)) {
        ASSERT_EQ(andnotbits[i], andnotbit);
    }
}

TEST(BitsetEnumerator, BitsetXor) {
    Bitset bitset1(200);
    defs::inds setbits1{1, 4, 7, 32, 89, 124, 199};
    Bitset bitset2(200);
    defs::inds setbits2{1, 4, 5, 7, 89, 91, 123};

    defs::inds xorbits{5, 32, 91, 123, 124, 199};

    bitset1.set(setbits1);
    bitset2.set(setbits2);

    BitsetXorEnumerator enumerator(bitset1, bitset2);

    size_t xorbit;
    size_t i = ~0ul;
    while (enumerator.next(xorbit, i)) {
        ASSERT_EQ(xorbits[i], xorbit);
    }
}

TEST(BitsetEnumerator, DeterminantSet) {
    Determinant det(50);
    defs::inds spinorbs{1, 4, 7, 32, 50};//, 54, 60, 89};
    det.set(spinorbs);

    OccupiedOrbitals occorbs(det);
    for (size_t i=0; i<occorbs.m_nind; ++i){
        std::cout << occorbs.m_inds[i] <<std::endl;
    }

}