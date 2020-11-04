//
// Created by RJA on 24/09/2020.
//

#include "gtest/gtest.h"
#include "src/core/basis/BosonOnvConnection.h"
#include "src/core/field/Elements.h"


TEST(BosonOnvConnection, NoChange) {
    size_t nmode = 4ul;

    elements::BosonOnv ket(nmode);
    elements::BosonOnv bra(nmode);

    ket = {2, 4, 0, 1};
    bra = {2, 4, 0, 1};

    BosonOnvConnection pc(ket, bra);

    ASSERT_EQ(pc.nchanged_mode(), 0);
}

TEST(BosonOnvConnection, SingleChange) {
    size_t nmode = 4ul;
    size_t occ_cutoff = 6ul;

    elements::BosonOnv ket(nmode);
    elements::BosonOnv bra(nmode);
    elements::BosonOnv work_bonv(nmode);

    for (size_t imode = 0; imode < nmode; ++imode) {
        for (size_t idelta = 1; idelta < occ_cutoff; ++idelta) {
            BosonOnvConnection pc(ket, bra);

            ket = {2, 4, 0, 1};
            bra = {2, 4, 0, 1};

            bra(imode) += idelta;
            pc.connect(bra, ket);
            ASSERT_EQ(pc.nchanged_mode(), 1);
            ASSERT_EQ(pc.changed_mode(0), imode);
            ASSERT_EQ(pc.changes(0), idelta);
            pc.apply(ket, work_bonv);
            ASSERT_EQ(work_bonv, bra);
        }
    }
}

TEST(BosonOnvConnection, DoubleChange) {
    size_t nmode = 4ul;
    size_t occ_cutoff = 6ul;

    elements::BosonOnv ket(nmode);
    elements::BosonOnv bra(nmode);
    elements::BosonOnv work_bonv(nmode);

    for (size_t imode1 = 0; imode1 < nmode; ++imode1) {
        for (size_t imode2 = imode1+1; imode2 < nmode; ++imode2) {
            for (size_t idelta1 = 1; idelta1 < occ_cutoff; ++idelta1) {
                for (size_t idelta2 = 1; idelta2 < occ_cutoff; ++idelta2) {

                    ket = {2, 4, 0, 1};
                    bra = {2, 4, 0, 1};

                    BosonOnvConnection pc(ket, bra);
                    bra(imode1) += idelta1;
                    bra(imode2) += idelta2;
                    pc.connect(bra, ket);

                    ASSERT_EQ(pc.nchanged_mode(), 2);
                    ASSERT_EQ(pc.changed_mode(0), imode1);
                    ASSERT_EQ(pc.changed_mode(1), imode2);
                    ASSERT_EQ(pc.changes(0), idelta1);
                    ASSERT_EQ(pc.changes(1), idelta2);
                    pc.apply(ket, work_bonv);
                    ASSERT_EQ(work_bonv, bra);
                }
            }
        }
    }
}
