//
// Created by rja on 13/06/22.
//

#include "gtest/gtest.h"
#include "M7_lib/util/Exsig.h"

TEST(UtilExsig, EncodeDecode) {
    using namespace exsig;
    /*
     * assert that all excitation signatures are encoded and decoded correctly
     */
    size_t exsig;
    for (size_t ncref = 0ul; ncref <= defs::exsig_nop_mask_frm; ++ncref) {
        for (size_t nannf = 0ul; nannf <= defs::exsig_nop_mask_frm; ++nannf) {
            for (size_t ncreb = 0ul; ncreb <= defs::exsig_nop_mask_bos; ++ncreb) {
                for (size_t nannb = 0ul; nannb <= defs::exsig_nop_mask_bos; ++nannb) {
                    exsig = encode(ncref, nannf, ncreb, nannb);
                    ASSERT_EQ(decode_nfrm_cre(exsig), ncref);
                    ASSERT_EQ(decode_nfrm_ann(exsig), nannf);
                    ASSERT_EQ(decode_nbos_cre(exsig), ncreb);
                    ASSERT_EQ(decode_nbos_ann(exsig), nannb);
                }
                exsig = encode(ncref, nannf, ncreb, defs::exsig_nop_mask_bos + 1);
                ASSERT_EQ(exsig, ~0ul);
            }
            exsig = encode(ncref, nannf, defs::exsig_nop_mask_bos + 1, 0ul);
            ASSERT_EQ(exsig, ~0ul);
        }
        exsig = encode(ncref, defs::exsig_nop_mask_frm + 1, 0ul, 0ul);
        ASSERT_EQ(exsig, ~0ul);
    }
    exsig = encode(defs::exsig_nop_mask_frm + 1, 0ul, 0ul, 0ul);
    ASSERT_EQ(exsig, ~0ul);
}

TEST(UtilExsig, RanksigContribs) {
    using namespace exsig;
    auto ranksig = encode(4, 7, 3, 2);
    ASSERT_NE(ranksig, ~0ul);
    ASSERT_EQ(ncontrib_frm(ranksig), 5); // 4732, 3632, 2532, 1432, 0332
    ASSERT_EQ(ncontrib_bos(ranksig), 3); // 4732, 4721, 4710
    ASSERT_TRUE(contribs_to(encode(4,7,3,2), ranksig));
    ASSERT_TRUE(contribs_to(encode(3,6,3,2), ranksig));
    ASSERT_TRUE(contribs_to(encode(2,5,3,2), ranksig));
    ASSERT_TRUE(contribs_to(encode(1,4,3,2), ranksig));
    ASSERT_TRUE(contribs_to(encode(0,3,3,2), ranksig));

    ASSERT_TRUE(contribs_to(encode(4,7,2,1), ranksig));
    ASSERT_TRUE(contribs_to(encode(3,6,2,1), ranksig));
    ASSERT_TRUE(contribs_to(encode(2,5,2,1), ranksig));
    ASSERT_TRUE(contribs_to(encode(1,4,2,1), ranksig));
    ASSERT_TRUE(contribs_to(encode(0,3,2,1), ranksig));

    ASSERT_TRUE(contribs_to(encode(4,7,1,0), ranksig));
    ASSERT_TRUE(contribs_to(encode(3,6,1,0), ranksig));
    ASSERT_TRUE(contribs_to(encode(2,5,1,0), ranksig));
    ASSERT_TRUE(contribs_to(encode(1,4,1,0), ranksig));
    ASSERT_TRUE(contribs_to(encode(0,3,1,0), ranksig));

    ASSERT_FALSE(contribs_to(0ul, ranksig));
}
#if 0
TEST(Utils, SetAllExsigsFromRanksig) {
    size_t ranksig;
    std::array<bool, defs::nexsig> exsigs{};

    ranksig = conn_utils::encode_exsig(4, 4, 1, 1);
    exsigs.fill(false);
    ASSERT_EQ(std::count(exsigs.cbegin(), exsigs.cend(), true), 0);
    conn_utils::add_exsigs(ranksig, exsigs);
    ASSERT_TRUE(exsigs[conn_utils::encode_exsig(4, 4, 1, 1)]);
    ASSERT_TRUE(exsigs[conn_utils::encode_exsig(3, 3, 1, 1)]);
    ASSERT_TRUE(exsigs[conn_utils::encode_exsig(2, 2, 1, 1)]);
    ASSERT_TRUE(exsigs[conn_utils::encode_exsig(1, 1, 1, 1)]);
    ASSERT_TRUE(exsigs[conn_utils::encode_exsig(0, 0, 1, 1)]);
    ASSERT_TRUE(exsigs[conn_utils::encode_exsig(4, 4, 0, 0)]);
    ASSERT_TRUE(exsigs[conn_utils::encode_exsig(3, 3, 0, 0)]);
    ASSERT_TRUE(exsigs[conn_utils::encode_exsig(2, 2, 0, 0)]);
    ASSERT_TRUE(exsigs[conn_utils::encode_exsig(1, 1, 0, 0)]);
    ASSERT_TRUE(exsigs[conn_utils::encode_exsig(0, 0, 0, 0)]);
    ASSERT_EQ(std::count(exsigs.cbegin(), exsigs.cend(), true), 10);

    ranksig = conn_utils::encode_exsig(4, 4, 0, 0);
    exsigs.fill(false);
    ASSERT_EQ(std::count(exsigs.cbegin(), exsigs.cend(), true), 0);
    conn_utils::add_exsigs(ranksig, exsigs);
    ASSERT_TRUE(exsigs[conn_utils::encode_exsig(4, 4, 0, 0)]);
    ASSERT_TRUE(exsigs[conn_utils::encode_exsig(3, 3, 0, 0)]);
    ASSERT_TRUE(exsigs[conn_utils::encode_exsig(2, 2, 0, 0)]);
    ASSERT_TRUE(exsigs[conn_utils::encode_exsig(1, 1, 0, 0)]);
    ASSERT_TRUE(exsigs[conn_utils::encode_exsig(0, 0, 0, 0)]);
    ASSERT_EQ(std::count(exsigs.cbegin(), exsigs.cend(), true), 5);

    ranksig = conn_utils::encode_exsig(3, 3, 1, 0);
    exsigs.fill(false);
    ASSERT_EQ(std::count(exsigs.cbegin(), exsigs.cend(), true), 0);
    conn_utils::add_exsigs(ranksig, exsigs);
    ASSERT_TRUE(exsigs[conn_utils::encode_exsig(3, 3, 1, 0)]);
    ASSERT_TRUE(exsigs[conn_utils::encode_exsig(2, 2, 1, 0)]);
    ASSERT_TRUE(exsigs[conn_utils::encode_exsig(1, 1, 1, 0)]);
    ASSERT_TRUE(exsigs[conn_utils::encode_exsig(0, 0, 1, 0)]);
    ASSERT_EQ(std::count(exsigs.cbegin(), exsigs.cend(), true), 4);


    ranksig = conn_utils::encode_exsig(2, 2, 0, 1);
    exsigs.fill(false);
    ASSERT_EQ(std::count(exsigs.cbegin(), exsigs.cend(), true), 0);
    conn_utils::add_exsigs(ranksig, exsigs);
    ASSERT_TRUE(exsigs[conn_utils::encode_exsig(2, 2, 0, 1)]);
    ASSERT_TRUE(exsigs[conn_utils::encode_exsig(1, 1, 0, 1)]);
    ASSERT_TRUE(exsigs[conn_utils::encode_exsig(0, 0, 0, 1)]);
    ASSERT_EQ(std::count(exsigs.cbegin(), exsigs.cend(), true), 3);
}
#endif