//
// Created by RJA on 24/09/2020.
//

#include <src/core/foreach/ForeachVirtual.h>
#include "gtest/gtest.h"
#include "src/core/table/BufferedFields.h"
#include "src/core/connection/Connections.h"

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
    size_t chk_occ_fac_square(const field::BosOnv& src, const defs::inds& icres, const defs::inds& ianns){
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

    bool one_conn(const field::BosOnv& src, defs::inds icres, defs::inds ianns){
        auto chk = chk_occ_fac_square(src, icres, ianns);
        conn::BosOnv conn(src.m_nmode);
        BosOps com(src.m_nmode);
        /*
         * remove common indices from the creation and annihilation vectors and place them in the common operators object
         */
        for (size_t imode=0ul; imode<src.m_nmode; ++imode){
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
            if (noccur) com.add({imode, noccur});
        }
        conn.m_cre.from_vector(icres);
        conn.m_ann.from_vector(ianns);
        DEBUG_ASSERT_EQ(conn.m_cre.to_vector(), icres, "incorrect connection encoding");
        DEBUG_ASSERT_EQ(conn.m_ann.to_vector(), ianns, "incorrect connection encoding");

        auto occ_fac = conn.occ_fac_square(src, com);
        //DEBUG_ASSERT_EQ(chk, occ_fac, "occ fac mismatch");
        return chk==occ_fac;
    }

    struct ForeachBase : foreach_virtual::rtnd::Ordered<false, true>{
        ForeachBase(size_t nmode, size_t nop): foreach_virtual::rtnd::Ordered<false, true>(nmode, nop){}
    };

    struct AnnForeach : ForeachBase {
        const field::BosOnv& m_src;
        const defs::inds& m_cre_inds;
        AnnForeach(size_t nmode, size_t nop, const defs::inds& cre_inds, const field::BosOnv& src):
            ForeachBase(nmode, nop), m_src(src), m_cre_inds(cre_inds){}

        void body(const foreach_virtual::rtnd::inds_t &value, size_t iiter) override {
            ASSERT_TRUE(one_conn(m_src, m_cre_inds, value));
        }
    };

    struct CreForeach : ForeachBase {
        const field::BosOnv& m_src;
        const size_t m_nmode, m_nann;
        CreForeach(size_t nmode, size_t ncre, size_t nann, const field::BosOnv& src):
            ForeachBase(nmode, ncre), m_src(src), m_nmode(nmode), m_nann(nann){}

        void body(const foreach_virtual::rtnd::inds_t &value, size_t iiter) override {
            AnnForeach foreach(m_nmode, m_nann, value, m_src);
            foreach.loop();
        }
    };
}

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

TEST(BosonOnvConnection, OccFac_0001) {
    const size_t nmode = 8;
    buffered::BosOnv mbf(nmode);
    mbf = {3, 4, 1, 2, 0, 5, 0, 2};
    conn::BosOnv conn(nmode);
    using namespace bos_onv_connection_test;
    CreForeach foreach(nmode, 0, 1, mbf);
    foreach.loop();
}

TEST(BosonOnvConnection, OccFac_0010) {
    const size_t nmode = 8;
    buffered::BosOnv mbf(nmode);
    mbf = {3, 4, 1, 2, 0, 5, 0, 2};
    conn::BosOnv conn(nmode);
    using namespace bos_onv_connection_test;
    CreForeach foreach(nmode, 1, 0, mbf);
    foreach.loop();
}
TEST(BosonOnvConnection, OccFac_0020) {
    const size_t nmode = 8;
    buffered::BosOnv mbf(nmode);
    mbf = {3, 4, 1, 2, 0, 5, 0, 2};
    conn::BosOnv conn(nmode);
    using namespace bos_onv_connection_test;
    CreForeach foreach(nmode, 2, 0, mbf);
    foreach.loop();
}
TEST(BosonOnvConnection, OccFac_0002) {
    const size_t nmode = 8;
    buffered::BosOnv mbf(nmode);
    mbf = {3, 4, 1, 2, 0, 5, 0, 2};
    conn::BosOnv conn(nmode);
    using namespace bos_onv_connection_test;
    CreForeach foreach(nmode, 0, 2, mbf);
    foreach.loop();
}
TEST(BosonOnvConnection, OccFac_0022) {
    const size_t nmode = 6;
    buffered::BosOnv mbf(nmode);
    mbf = {3, 4, 1, 2, 0, 5};
    conn::BosOnv conn(nmode);
    using namespace bos_onv_connection_test;
    CreForeach foreach(nmode, 2, 2, mbf);
    foreach.loop();
}
TEST(BosonOnvConnection, OccFac_0012) {
    const size_t nmode = 6;
    buffered::BosOnv mbf(nmode);
    mbf = {3, 4, 1, 2, 0, 5};
    conn::BosOnv conn(nmode);
    using namespace bos_onv_connection_test;
    CreForeach foreach(nmode, 1, 2, mbf);
    foreach.loop();
}
TEST(BosonOnvConnection, OccFac_0033) {
    const size_t nmode = 5;
    buffered::BosOnv mbf(nmode);
    mbf = {3, 4, 0, 2, 5};
    conn::BosOnv conn(nmode);
    using namespace bos_onv_connection_test;
    CreForeach foreach(nmode, 3, 3, mbf);
    foreach.loop();
}
TEST(BosonOnvConnection, OccFac_0032) {
    const size_t nmode = 5;
    buffered::BosOnv mbf(nmode);
    mbf = {3, 4, 0, 2, 5};
    conn::BosOnv conn(nmode);
    using namespace bos_onv_connection_test;
    CreForeach foreach(nmode, 3, 2, mbf);
    foreach.loop();
}
TEST(BosonOnvConnection, OccFac_0031) {
    const size_t nmode = 5;
    buffered::BosOnv mbf(nmode);
    mbf = {3, 4, 0, 2, 5};
    conn::BosOnv conn(nmode);
    using namespace bos_onv_connection_test;
    CreForeach foreach(nmode, 3, 1, mbf);
    foreach.loop();
}