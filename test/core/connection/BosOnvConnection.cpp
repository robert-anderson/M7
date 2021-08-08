//
// Created by RJA on 24/09/2020.
//

#include "gtest/gtest.h"
#include "src/core/table/BufferedFields.h"
#include "src/core/connection/Connections.h"


TEST(BosonOnvConnection, NoChange) {
    size_t nmode = 4ul;

    buffered::BosOnv ket(nmode);
    buffered::BosOnv bra(nmode);

    ket = defs::inds{2, 4, 0, 1};
    bra =  defs::inds{2, 4, 0, 1};

    conn::BosOnv pc(nmode);
    pc.connect(ket, bra);

    ASSERT_EQ(pc.size(), 0);
}

TEST(BosonOnvConnection, SingleChange) {
    size_t nmode = 4ul;
    size_t occ_cutoff = 6ul;

    buffered::BosOnv ket(nmode);
    buffered::BosOnv bra(nmode);
    buffered::BosOnv work_bonv(nmode);

    for (size_t imode = 0; imode < nmode; ++imode) {
        for (size_t idelta = 1; idelta < occ_cutoff; ++idelta) {
            BosOnvConnection conn(nmode);

            ket = {2, 4, 0, 1};
            bra = {2, 4, 0, 1};

            bra[imode] += idelta;
            conn.connect(ket, bra);
            ASSERT_EQ(conn.m_cre.size(), idelta);
            ASSERT_EQ(conn.m_ann.size(), 0);
            ASSERT_EQ(conn.size(), idelta);
            ASSERT_EQ(conn.m_cre[0].m_imode, imode);
            ASSERT_EQ(conn.m_cre[0].m_nop, idelta);
            conn.apply(ket, work_bonv);
            ASSERT_EQ(work_bonv, bra);
        }
    }
}

TEST(BosonOnvConnection, DoubleChange) {
    size_t nmode = 4ul;
    size_t occ_cutoff = 6ul;

    buffered::BosOnv ket(nmode);
    buffered::BosOnv bra(nmode);
    buffered::BosOnv work_bonv(nmode);

    for (size_t imode0 = 0; imode0 < nmode; ++imode0) {
        for (size_t imode1 = imode0 + 1; imode1 < nmode; ++imode1) {
            for (size_t idelta0 = 1; idelta0 < occ_cutoff; ++idelta0) {
                for (size_t idelta1 = 1; idelta1 < occ_cutoff; ++idelta1) {

                    ket = {2, 4, 0, 1};
                    bra = {2, 4, 0, 1};

                    BosOnvConnection pc(nmode);
                    bra[imode0] += idelta0;
                    bra[imode1] += idelta1;
                    pc.connect(ket, bra);

                    ASSERT_EQ(pc.size(), idelta0+idelta1);
                    ASSERT_EQ(pc.m_cre[0].m_imode, imode0);
                    ASSERT_EQ(pc.m_cre[1].m_imode, imode1);
                    ASSERT_EQ(pc.m_cre[0].m_nop, idelta0);
                    ASSERT_EQ(pc.m_cre[1].m_nop, idelta1);
                    pc.apply(ket, work_bonv);
                    ASSERT_EQ(work_bonv, bra);
                }
            }
        }
    }
}
