//
// Created by RJA on 24/09/2020.
//

#include "gtest/gtest.h"
#include "src/core/basis/PermanentConnection.h"
#include "src/core/basis/Permanent.h"

// todo setup function



TEST(PermanentConnection, NoChange) {
    size_t nmode = 4ul;
    size_t occ_cutoff = 6ul;

    Permanent ket(nmode, occ_cutoff);
    Permanent bra(nmode, occ_cutoff);

    ket = {2,4,0,1};
    bra = {2,4,0,1};

    PermanentConnection pc(ket, bra);

    ASSERT_EQ(pc.nchanged_mode(), 0);
}

TEST(PermanentConnection, SingleChange){
    size_t nmode = 4ul;
    size_t occ_cutoff = 6ul;

    Permanent ket(nmode, occ_cutoff);
    Permanent bra(nmode, occ_cutoff);

    for(size_t imode = 0; imode < nmode; ++imode){
        for (size_t idelta = 1; idelta < occ_cutoff; ++idelta) {
            PermanentConnection pc(ket, bra);

            ket = {2,4,0,1};
            bra = {2,4,0,1};

            bra(imode) += idelta;
            pc.connect(bra, ket);

            ASSERT_EQ(pc.nchanged_mode(), 1);
            ASSERT_EQ(pc.changed_modes(0), imode);
            ASSERT_EQ(pc.changes(0), idelta);
            }
        }
}


TEST(PermanentConnection, DoubleChange){
    size_t nmode = 4ul;
    size_t occ_cutoff = 6ul;

    Permanent ket(nmode, occ_cutoff);
    Permanent bra(nmode, occ_cutoff);

    ket = {2,4,0,1};
    bra = {2,4,0,1};

    //std::cout << "testing initial boson occupation vector with elements:" << std::endl;
    for (auto& item : ket.to_vector()) std::cout << (int)item << std::endl;

    for(size_t imode1 = 0; imode1 < nmode; ++imode1){
        for(size_t imode2 = 0; imode2 < imode1; ++imode2){
            for (size_t idelta1 = 1; idelta1 < occ_cutoff; ++idelta1) {
                for (size_t idelta2 = 1; idelta2 < occ_cutoff; ++idelta2) {
                    PermanentConnection pc(ket, bra);
                    bra(imode1) += idelta1;
                    bra(imode2) += idelta2;
                    pc.connect(bra, ket);

                    ASSERT_EQ(pc.nchanged_mode(), 2);
                    ASSERT_EQ(pc.changed_modes(0), imode1);
                    ASSERT_EQ(pc.changed_modes(1), imode2);
                    ASSERT_EQ(pc.changes(0), idelta1);
                    ASSERT_EQ(pc.changes(1), idelta2);
                }
        }
    }
}
}