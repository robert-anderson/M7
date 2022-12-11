//
// Created by rja on 13/06/22.
//

#include <algorithm>
#include "gtest/gtest.h"
#include "M7_lib/connection/OpSig.h"

TEST(UtilExsig, EncodeDecode) {
    /*
     * assert that all excitation signatures are encoded and decoded correctly
     */
    for (uint_t ncref = 0ul; ncref <= opsig::c_nop_mask_frm; ++ncref) {
        for (uint_t nannf = 0ul; nannf <= opsig::c_nop_mask_frm; ++nannf) {
            for (uint_t ncreb = 0ul; ncreb <= opsig::c_nop_mask_bos; ++ncreb) {
                for (uint_t nannb = 0ul; nannb <= opsig::c_nop_mask_bos; ++nannb) {
                    OpSig opsig({ncref, nannf}, {ncreb, nannb});
                    ASSERT_EQ(opsig.nfrm_cre(), ncref);
                    ASSERT_EQ(opsig.nfrm_ann(), nannf);
                    ASSERT_EQ(opsig.nbos_cre(), ncreb);
                    ASSERT_EQ(opsig.nbos_ann(), nannb);
                }
                OpSig opsig({ncref, nannf}, {ncreb, opsig::c_nop_mask_bos + 1});
                ASSERT_FALSE(opsig.is_valid());
            }
            OpSig opsig({ncref, nannf}, {opsig::c_nop_mask_bos + 1, 0ul});
            ASSERT_FALSE(opsig.is_valid());
        }
        OpSig opsig({ncref, opsig::c_nop_mask_frm + 1}, {0ul, 0ul});
        ASSERT_FALSE(opsig.is_valid());
    }
    OpSig opsig({opsig::c_nop_mask_frm + 1, 0ul}, {0ul, 0ul});
    ASSERT_FALSE(opsig.is_valid());
}

TEST(UtilExsig, RanksigContribs) {
    OpSig ranksig({4ul, 7ul}, {3ul, 2ul});
    ASSERT_TRUE(ranksig.is_valid());
    ASSERT_EQ(ranksig.ncontrib_frm(), 5ul); // 4732, 3632, 2532, 1432, 0332
    ASSERT_EQ(ranksig.ncontrib_bos(), 3ul); // 4732, 4721, 4710

    ASSERT_TRUE(ranksig.takes_contribs_from({{4, 7}, {3, 2}}));
    ASSERT_TRUE(ranksig.takes_contribs_from({{3, 6}, {3, 2}}));
    ASSERT_TRUE(ranksig.takes_contribs_from({{2, 5}, {3, 2}}));
    ASSERT_TRUE(ranksig.takes_contribs_from({{1, 4}, {3, 2}}));
    ASSERT_TRUE(ranksig.takes_contribs_from({{0, 3}, {3, 2}}));

    ASSERT_TRUE(ranksig.takes_contribs_from({{4, 7}, {2, 1}}));
    ASSERT_TRUE(ranksig.takes_contribs_from({{3, 6}, {2, 1}}));
    ASSERT_TRUE(ranksig.takes_contribs_from({{2, 5}, {2, 1}}));
    ASSERT_TRUE(ranksig.takes_contribs_from({{1, 4}, {2, 1}}));
    ASSERT_TRUE(ranksig.takes_contribs_from({{0, 3}, {2, 1}}));

    ASSERT_TRUE(ranksig.takes_contribs_from({{4, 7}, {1, 0}}));
    ASSERT_TRUE(ranksig.takes_contribs_from({{3, 6}, {1, 0}}));
    ASSERT_TRUE(ranksig.takes_contribs_from({{2, 5}, {1, 0}}));
    ASSERT_TRUE(ranksig.takes_contribs_from({{1, 4}, {1, 0}}));
    ASSERT_TRUE(ranksig.takes_contribs_from({{0, 3}, {1, 0}}));

    ASSERT_FALSE(ranksig.takes_contribs_from(opsig::c_zero));
}