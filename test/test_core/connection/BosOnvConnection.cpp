//
// Created by Robert J. Anderson on 24/09/2020.
//

#include "gtest/gtest.h"
#include "M7_lib/table/BufferedFields.h"
#include "M7_lib/connection/Connections.h"

namespace bos_onv_connection_test {
    /**
     * apply the creation and annihilation operators (with common mode indices included) to the source BosOnv by
     * modifying a working BosOnv each time an operator is applied. the operator index vectors differ from those found
     * in a boson ONV connection in that repetitions are not grouped by mode index, and creation and annihilation
     * operators are allowed to have indices in common.
     *
     * this is a slow sanity check on the fast BosOnvConnection::occ_fac_square method.
     *
     * @param src
     *  boson ONV to be acted on with the normal ordered SQ operator product
     * @param icres
     *  creation operator product
     * @param ianns
     *  annihilation operator product
     * @return
     *  square of the occupation factor
     */
    size_t chk_occ_fac_square(const field::BosOnv& src, const defs::ivec_t& icres, const defs::ivec_t& ianns){
        buffered::BosOnv work(src);
        size_t tot = 1;
        for (auto& iann : ianns) {
            tot*=work[iann];
            if (work[iann]) work[iann]--;
        }
        for (auto& icre : icres) {
            tot*=work[icre]+1;
            work[icre]++;
        }
        return tot;
    }

    bool one_conn(const field::BosOnv& src, defs::ivec_t icres, defs::ivec_t ianns){
        auto chk = chk_occ_fac_square(src, icres, ianns);
        conn::BosOnv conn(src.m_basis.m_nmode);
        BosOps com(src.m_basis.m_nmode);
        /*
         * remove common indices from the creation and annihilation vectors and place them in the common operators object
         */
        for (size_t imode=0ul; imode<src.m_basis.m_nmode; ++imode){
            size_t noccur = 0;
            while (true) {
                auto cre_it = std::find(icres.begin(), icres.end(), imode);
                if (cre_it==icres.end()) break;
                auto ann_it = std::find(ianns.begin(), ianns.end(), imode);
                if (ann_it==ianns.end()) break;
                icres.erase(cre_it);
                ianns.erase(ann_it);
                ++noccur;
            }
            if (noccur) com.add(imode, noccur);
        }
        conn.m_cre.set(icres);
        conn.m_ann.set(ianns);
        DEBUG_ASSERT_EQ(conn.m_cre.get(), icres, "incorrect connection encoding");
        DEBUG_ASSERT_EQ(conn.m_ann.get(), ianns, "incorrect connection encoding");

        auto occ_fac = src.occ_fac_square(conn, com);
        if (!conn.size()) {
            // diagonal occ facs have another method which doesn't require the BosOnvConnection:
            return (chk==occ_fac) && (chk==src.occ_fac_square(com));
        }
        return chk==occ_fac;
    }

    void foreach(size_t nmode, size_t ncre, size_t nann, const field::BosOnv& src) {
        using namespace basic_foreach::rtnd;
        auto cre_foreach = basic_foreach::rtnd::Ordered<false, true>(nmode, ncre);
        auto fn = [&](const inds_t& cre_inds) {
            auto ann_foreach = basic_foreach::rtnd::Ordered<false, true>(nmode, nann);
            auto fn = [&](const inds_t& ann_inds) {
                ASSERT_TRUE(one_conn(src, cre_inds, ann_inds));
            };
            ann_foreach.loop(fn);            
        };
        cre_foreach.loop(fn);
    }
}

TEST(BosonOnvConnection, NoChange) {
    size_t nmode = 4ul;

    buffered::BosOnv ket(nmode);
    buffered::BosOnv bra(nmode);

    ket = defs::ivec_t{2, 4, 0, 1};
    bra =  defs::ivec_t{2, 4, 0, 1};

    conn::BosOnv pc(ket);
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
            BosOnvConnection conn(ket);

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

                    BosOnvConnection pc(ket);
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

TEST(BosonOnvConnection, OccFac_0000) {
    const size_t nmode = 8;
    buffered::BosOnv mbf(nmode);
    mbf = {3, 4, 1, 2, 0, 5, 0, 2};
    conn::BosOnv conn(nmode);
    bos_onv_connection_test::foreach(nmode, 0, 0, mbf);
}

TEST(BosonOnvConnection, OccFac_0001) {
    const size_t nmode = 8;
    buffered::BosOnv mbf(nmode);
    mbf = {3, 4, 1, 2, 0, 5, 0, 2};
    conn::BosOnv conn(mbf);
    bos_onv_connection_test::foreach(nmode, 0, 1, mbf);
}

TEST(BosonOnvConnection, OccFac_0010) {
    const size_t nmode = 8;
    buffered::BosOnv mbf(nmode);
    mbf = {3, 4, 1, 2, 0, 5, 0, 2};
    conn::BosOnv conn(mbf);
    bos_onv_connection_test::foreach(nmode, 1, 0, mbf);
}
TEST(BosonOnvConnection, OccFac_0020) {
    const size_t nmode = 8;
    buffered::BosOnv mbf(nmode);
    mbf = {3, 4, 1, 2, 0, 5, 0, 2};
    conn::BosOnv conn(mbf);
    bos_onv_connection_test::foreach(nmode, 2, 0, mbf);
}
TEST(BosonOnvConnection, OccFac_0002) {
    const size_t nmode = 8;
    buffered::BosOnv mbf(nmode);
    mbf = {3, 4, 1, 2, 0, 5, 0, 2};
    conn::BosOnv conn(mbf);
    bos_onv_connection_test::foreach(nmode, 0, 2, mbf);
}
TEST(BosonOnvConnection, OccFac_0022) {
    const size_t nmode = 6;
    buffered::BosOnv mbf(nmode);
    mbf = {3, 4, 1, 2, 0, 5};
    conn::BosOnv conn(mbf);
    bos_onv_connection_test::foreach(nmode, 2, 2, mbf);
}
TEST(BosonOnvConnection, OccFac_0012) {
    const size_t nmode = 6;
    buffered::BosOnv mbf(nmode);
    mbf = {3, 4, 1, 2, 0, 5};
    conn::BosOnv conn(mbf);
    bos_onv_connection_test::foreach(nmode, 1, 2, mbf);
}
TEST(BosonOnvConnection, OccFac_0033) {
    const size_t nmode = 5;
    buffered::BosOnv mbf(nmode);
    mbf = {3, 4, 0, 2, 5};
    conn::BosOnv conn(mbf);
    bos_onv_connection_test::foreach(nmode, 3, 3, mbf);
}
TEST(BosonOnvConnection, OccFac_0032) {
    const size_t nmode = 5;
    buffered::BosOnv mbf(nmode);
    mbf = {3, 4, 0, 2, 5};
    conn::BosOnv conn(mbf);
    bos_onv_connection_test::foreach(nmode, 3, 2, mbf);
}
TEST(BosonOnvConnection, OccFac_0031) {
    const size_t nmode = 5;
    buffered::BosOnv mbf(nmode);
    mbf = {3, 4, 0, 2, 5};
    conn::BosOnv conn(mbf);
    bos_onv_connection_test::foreach(nmode, 3, 1, mbf);
}
